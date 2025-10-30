#include <bits/stdc++.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <chrono>
#include <omp.h>  // OpenMP支持

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// 时间统计工具
using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;
inline double elapsed_ms(TimePoint start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

using Point = bg::model::d2::point_xy<int>;
using Polygon = bg::model::polygon<Point>;
using Box = bg::model::box<Point>;
using RTreeValue = std::pair<Box, int>;
using RTree = bgi::rtree<RTreeValue, bgi::rstar<24>>;

struct PhaseStopwatch {
    std::chrono::high_resolution_clock::time_point t0;
    double& sink;
    bool active;
    PhaseStopwatch(const char* name, double& s)
        : sink(s), active(name != nullptr) {
        if (active) {
            t0 = std::chrono::high_resolution_clock::now();
        }
    }
    ~PhaseStopwatch() {
        if (active) {
            auto t1 = std::chrono::high_resolution_clock::now();
            sink += std::chrono::duration<double>(t1 - t0).count();
        }
    }
};

double g_time_parse = 0.0;
double g_time_pre = 0.0;
double g_time_s1 = 0.0;
double g_time_s2 = 0.0;
double g_time_slice = 0.0;
double g_time_output = 0.0;
bool   g_profile = true;

// 多边形信息结构体
struct PolygonInfo {
    int id;
    int layer;
    std::vector<std::pair<int, int>> vertices;
    Polygon geom;
    Box bbox;
    std::vector<std::array<int, 3>> horizEdges;
    std::vector<std::array<int, 3>> vertEdges;
    int originalId;
};

// Gate边结构
struct GateEdge {
    int to;
    int polyId;
};

// ==================== 线程控制 ====================
int g_threadCount = 1;

inline int get_thread_count() {
    return g_threadCount;
}

inline void set_omp_threads(int n) {
    g_threadCount = std::max(1, n);
#ifdef _OPENMP
    omp_set_num_threads(g_threadCount);
#endif
}

// ==================== 全局变量 ====================
std::vector<PolygonInfo> polygons;
std::unordered_map<int, std::vector<int>> layerPolygons;
std::unordered_map<int, RTree> layerIndex;
std::vector<std::vector<int>> adjacency;
std::vector<std::vector<GateEdge>> gateAdjacency;
std::vector<char> visited;
std::map<std::string, int> layerNameToId;
std::map<int, std::string> layerIdToName;
int nextLayerId = 1;
static std::vector<std::vector<int>> g_layerPolyIds;

struct GridIndex {
    int cell = 0;
    int64_t xmin = 0;
    int64_t ymin = 0;
    int64_t xmax = 0;
    int64_t ymax = 0;
    std::unordered_map<uint64_t, std::vector<int>> buckets;
};

static std::vector<GridIndex> g_grid;
static std::vector<int> g_seen_stamp;
static int g_seen_tick = 1;

inline uint64_t cell_key(int ix, int iy) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(ix)) << 32) ^ static_cast<uint32_t>(iy);
}

inline void ensure_seen_stamp() {
    if (g_seen_stamp.size() != polygons.size()) {
        g_seen_stamp.assign(polygons.size(), 0);
        g_seen_tick = 1;
    }
}

inline void grid_put(GridIndex& G, const Box& b, int pid) {
    if (G.cell <= 0) return;
    long long xmin = b.min_corner().x();
    long long ymin = b.min_corner().y();
    long long xmax = b.max_corner().x();
    long long ymax = b.max_corner().y();
    long long dx0 = std::max<long long>(0, xmin - G.xmin);
    long long dy0 = std::max<long long>(0, ymin - G.ymin);
    long long dx1 = std::max<long long>(0, xmax - G.xmin);
    long long dy1 = std::max<long long>(0, ymax - G.ymin);
    int ix0 = static_cast<int>(dx0 / G.cell);
    int iy0 = static_cast<int>(dy0 / G.cell);
    int ix1 = static_cast<int>(dx1 / G.cell);
    int iy1 = static_cast<int>(dy1 / G.cell);
    for (int ix = ix0; ix <= ix1; ++ix) {
        for (int iy = iy0; iy <= iy1; ++iy) {
            G.buckets[cell_key(ix, iy)].push_back(pid);
        }
    }
}

inline int estimate_cell(const std::vector<int>& ids) {
    if (ids.empty()) return 1024;
    std::vector<int> ws;
    std::vector<int> hs;
    ws.reserve(ids.size());
    hs.reserve(ids.size());
    for (int pid : ids) {
        const auto& box = polygons[pid].bbox;
        int w = std::max(1, box.max_corner().x() - box.min_corner().x());
        int h = std::max(1, box.max_corner().y() - box.min_corner().y());
        ws.push_back(w);
        hs.push_back(h);
    }
    auto median = [](std::vector<int>& vals) {
        size_t mid = vals.size() / 2;
        std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
        return std::max(1, vals[mid]);
    };
    int w = median(ws);
    int h = median(hs);
    return std::max(1, (w + h) / 2);
}

inline void build_uniform_grid_indices() {
    g_grid.clear();
    g_grid.resize(g_layerPolyIds.size());
    for (int L = 0; L < static_cast<int>(g_layerPolyIds.size()); ++L) {
        const auto& ids = g_layerPolyIds[L];
        if (ids.empty()) continue;
        GridIndex G;
        G.xmin = G.ymin = std::numeric_limits<int64_t>::max();
        G.xmax = G.ymax = std::numeric_limits<int64_t>::min();
        for (int pid : ids) {
            const auto& box = polygons[pid].bbox;
            G.xmin = std::min<int64_t>(G.xmin, box.min_corner().x());
            G.ymin = std::min<int64_t>(G.ymin, box.min_corner().y());
            G.xmax = std::max<int64_t>(G.xmax, box.max_corner().x());
            G.ymax = std::max<int64_t>(G.ymax, box.max_corner().y());
        }
        G.cell = estimate_cell(ids);
        G.buckets.reserve(ids.size() * 2);
        for (int pid : ids) {
            grid_put(G, polygons[pid].bbox, pid);
        }
        g_grid[L] = std::move(G);
    }
    ensure_seen_stamp();
}

struct OutPoly {
    int layer;
    std::vector<std::pair<int, int>> verts;
};

struct FastWriter {
    std::vector<char> buf;
    explicit FastWriter(size_t cap = 4u << 20) { buf.reserve(cap); }
    inline void put(char c) { buf.push_back(c); }
    inline void puts(const std::string& s) { buf.insert(buf.end(), s.begin(), s.end()); }
    inline void puti(int x) {
        long long val = x;
        if (val == 0) { buf.push_back('0'); return; }
        if (val < 0) { buf.push_back('-'); val = -val; }
        char tmp[32]; int n = 0;
        while (val > 0) { tmp[n++] = char('0' + (val % 10)); val /= 10; }
        while (n--) buf.push_back(tmp[n]);
    }
    inline void sp() { buf.push_back(' '); }
    inline void nl() { buf.push_back('\n'); }
    bool write_to_file(const std::string& path) {
        FILE* f = fopen(path.c_str(), "wb");
        if (!f) return false;
        size_t written = fwrite(buf.data(), 1, buf.size(), f);
        fclose(f);
        return written == buf.size();
    }
};

// 相交判断缓存（仅用于问题1优化）
// 注意：问题2和3构建全图时会产生大量缓存条目，反而降低性能
// 因此只在问题1的 namespace 内局部启用

// 命令行参数
std::string layoutFile, ruleFile, outputFile;
int threadCount = 1;
bool g_fastgeo = false;
bool g_s1bfs = true;
bool g_s2bfs = true;
int  g_slice_threads = 8;
bool g_slice_enable = true;

// 解析层名转换为层ID
int getLayerId(const std::string& layerName) {
    if (layerNameToId.find(layerName) == layerNameToId.end()) {
        layerNameToId[layerName] = nextLayerId;
        layerIdToName[nextLayerId] = layerName;
        nextLayerId++;
    }
    return layerNameToId[layerName];
}

// 解析命令行参数
bool parseCommandLine(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-layout" && i + 1 < argc) {
            layoutFile = argv[++i];
        } else if (arg == "-rule" && i + 1 < argc) {
            ruleFile = argv[++i];
        } else if (arg == "-output" && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (arg == "-thread" && i + 1 < argc) {
            threadCount = std::atoi(argv[++i]);
            if (threadCount < 1 || threadCount > 64) {
                std::cerr << "Error: thread count must be between 1 and 64\n";
                return false;
            }
        } else if (arg == "-fastgeo" && i + 1 < argc) {
            int v = std::atoi(argv[++i]);
            g_fastgeo = (v != 0);
        } else if (arg == "-s1bfs" && i + 1 < argc) {
            int v = std::atoi(argv[++i]);
            g_s1bfs = (v != 0);
        } else if (arg == "-s2bfs" && i + 1 < argc) {
            int v = std::atoi(argv[++i]);
            g_s2bfs = (v != 0);
        } else if (arg == "-slicethreads" && i + 1 < argc) {
            g_slice_threads = std::max(1, std::atoi(argv[++i]));
        } else if (arg == "-slice" && i + 1 < argc) {
            int v = std::atoi(argv[++i]);
            g_slice_enable = (v != 0);
        } else if (arg == "-profile" && i + 1 < argc) {
            g_profile = (std::atoi(argv[++i]) != 0);
        }
    }

    if (layoutFile.empty() || ruleFile.empty() || outputFile.empty()) {
        std::cerr << "Usage: " << argv[0]
                  << " -layout <file> -rule <file> -output <file>"
                  << " [-thread <n>] [-fastgeo 0|1] [-s1bfs 0|1] [-s2bfs 0|1]"
                  << " [-slice 0|1] [-slicethreads <n>] [-profile 0|1]\n";
        return false;
    }
    
    set_omp_threads(threadCount);
    
    std::cerr << "Configuration:\n";
    std::cerr << "  Layout: " << layoutFile << "\n";
    std::cerr << "  Rule: " << ruleFile << "\n";
    std::cerr << "  Output: " << outputFile << "\n";
    std::cerr << "  Threads: " << g_threadCount << "\n";
    std::cerr << "  Fast-geo: " << (g_fastgeo ? "ON" : "OFF") << "\n";
    std::cerr << "  S1 OnDemand BFS: " << (g_s1bfs ? "ON" : "OFF") << "\n";
    std::cerr << "  S2 OnDemand BFS: " << (g_s2bfs ? "ON" : "OFF") << "\n";
    std::cerr << "  Slice: " << (g_slice_enable ? "ON" : "OFF")
              << ", threads=" << g_slice_threads << "\n";
    std::cerr << "  Profile: " << (g_profile ? "ON" : "OFF") << "\n";

    return true;
}

// 手动解析整数
inline int fast_atoi(const char* str, const char** endptr) {
    int val = 0;
    bool negative = false;
    
    while (*str == ' ' || *str == '\t') ++str;
    
    if (*str == '-') {
        negative = true;
        ++str;
    } else if (*str == '+') {
        ++str;
    }
    
    while (*str >= '0' && *str <= '9') {
        val = val * 10 + (*str - '0');
        ++str;
    }
    
    if (endptr) *endptr = str;
    return negative ? -val : val;
}

inline void precompute_edges(PolygonInfo& poly) {
    poly.horizEdges.clear();
    poly.vertEdges.clear();
    const auto& V = poly.vertices;
    if (V.empty()) return;

    auto add_h = [&](int y, int x1, int x2) {
        if (x1 > x2) std::swap(x1, x2);
        poly.horizEdges.push_back({y, x1, x2});
    };
    auto add_v = [&](int x, int y1, int y2) {
        if (y1 > y2) std::swap(y1, y2);
        poly.vertEdges.push_back({x, y1, y2});
    };

    for (size_t i = 0; i < V.size(); ++i) {
        const auto& a = V[i];
        const auto& b = V[(i + 1) % V.size()];
        if (a.second == b.second) {
            add_h(a.second, a.first, b.first);
        } else if (a.first == b.first) {
            add_v(a.first, a.second, b.second);
        } else {
            // 非轴对齐边在题目设定下不应出现
        }
    }

    std::sort(poly.horizEdges.begin(), poly.horizEdges.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs[0] != rhs[0]) return lhs[0] < rhs[0];
        return lhs[1] < rhs[1];
    });
    std::sort(poly.vertEdges.begin(), poly.vertEdges.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs[0] != rhs[0]) return lhs[0] < rhs[0];
        return lhs[1] < rhs[1];
    });
}

// 解析版图文件
bool parse_layout_file(const std::string& filename) {
    std::ios_base::sync_with_stdio(false);
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open layout file: " << filename << "\n";
        return false;
    }
    
    polygons.reserve(1000000);
    
    std::string line;
    line.reserve(4096);
    int currentLayer = -1;
    int polygonId = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        const char* ptr = line.c_str();
        while (*ptr == ' ' || *ptr == '\t' || *ptr == '\r') ++ptr;
        if (*ptr == '\0') continue;
        
        if (*ptr == 'L' && (ptr[1] >= '0' && ptr[1] <= '9')) {
            bool hasParenthesis = false;
            for (const char* p = ptr; *p; ++p) {
                if (*p == '(') {
                    hasParenthesis = true;
                    break;
                }
            }
            
            if (!hasParenthesis) {
                currentLayer = getLayerId(ptr);
                continue;
            }
        }
        
        if (currentLayer != -1 && *ptr == '(') {
            std::vector<std::pair<int, int>> vertices;
            vertices.reserve(32);
            
            int minX = INT_MAX, maxX = INT_MIN;
            int minY = INT_MAX, maxY = INT_MIN;
            
            while (*ptr) {
                while (*ptr && *ptr != '(') ++ptr;
                if (!*ptr) break;
                ++ptr;
                
                const char* endptr;
                int x = fast_atoi(ptr, &endptr);
                ptr = endptr;
                
                while (*ptr && *ptr != ',') ++ptr;
                if (!*ptr) break;
                ++ptr;
                
                int y = fast_atoi(ptr, &endptr);
                ptr = endptr;
                
                while (*ptr && *ptr != ')') ++ptr;
                if (*ptr) ++ptr;
                
                vertices.push_back({x, y});
                
                if (x < minX) minX = x;
                if (x > maxX) maxX = x;
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
                
                while (*ptr && *ptr != '(' && *ptr != '\0') ++ptr;
            }
            
            if (!vertices.empty()) {
                PolygonInfo poly;
                poly.id = polygonId++;
                poly.layer = currentLayer;
                poly.vertices = std::move(vertices);
                poly.bbox = Box(Point(minX, minY), Point(maxX, maxY));
                
                for (const auto& v : poly.vertices) {
                    bg::append(poly.geom.outer(), Point(v.first, v.second));
                }
                
                if (!poly.vertices.empty() &&
                    poly.vertices.front() != poly.vertices.back()) {
                    const auto& first = poly.vertices[0];
                    bg::append(poly.geom.outer(), Point(first.first, first.second));
                }

                bg::correct(poly.geom);
                precompute_edges(poly);

                polygons.push_back(std::move(poly));
                layerPolygons[currentLayer].push_back(polygonId - 1);
            }
        }
    }
    
    file.close();
    
    std::cerr << "Building R-tree indexes (bulk loading)...\n";
    for (const auto& layerEntry : layerPolygons) {
        int layer = layerEntry.first;
        std::vector<RTreeValue> values;
        values.reserve(layerEntry.second.size());
        
        for (int id : layerEntry.second) {
            values.push_back(std::make_pair(polygons[id].bbox, id));
        }
        
        layerIndex[layer] = RTree(values.begin(), values.end());
    }
    std::cerr << "R-tree indexes built.\n";
    
    return true;
}

// 解析坐标
std::pair<int, int> parseCoordinate(const std::string& coordStr) {
    std::string clean = coordStr;
    clean.erase(std::remove(clean.begin(), clean.end(), '('), clean.end());
    clean.erase(std::remove(clean.begin(), clean.end(), ')'), clean.end());
    
    size_t commaPos = clean.find(',');
    if (commaPos != std::string::npos) {
        int x = std::stoi(clean.substr(0, commaPos));
        int y = std::stoi(clean.substr(commaPos + 1));
        return {x, y};
    }
    return {0, 0};
}

// 解析规则文件
bool parse_rule_file(const std::string& filename,
                     std::vector<std::pair<int, std::pair<int, int>>>& startPoints,
                     std::vector<std::vector<int>>& viaRules,
                     int& polyLayer,
                     int& aaLayer) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open rule file: " << filename << "\n";
        return false;
    }
    
    std::string line;
    bool inStartPos = false;
    bool inVia = false;
    bool inGate = false;
    
    polyLayer = -1;
    aaLayer = -1;
    
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty()) continue;
        
        if (line == "StartPos") {
            inStartPos = true;
            inVia = false;
            inGate = false;
            continue;
        }
        
        if (line == "Via") {
            inStartPos = false;
            inVia = true;
            inGate = false;
            continue;
        }
        
        if (line == "Gate") {
            inStartPos = false;
            inVia = false;
            inGate = true;
            continue;
        }
        
        if (inVia) {
            std::istringstream iss(line);
            std::vector<int> layers;
            std::string layerName;
            while (iss >> layerName) {
                layers.push_back(getLayerId(layerName));
            }
            
            if (!layers.empty()) {
                viaRules.push_back(layers);
            }
            continue;
        }
        
        if (inGate) {
            std::istringstream iss(line);
            std::string layer1, layer2;
            iss >> layer1 >> layer2;
            
            polyLayer = getLayerId(layer1);
            aaLayer = getLayerId(layer2);
            inGate = false;
            continue;
        }
        
        if (inStartPos) {
            std::istringstream iss(line);
            std::string layerName, coordStr;
            iss >> layerName >> coordStr;
            
            int layer = getLayerId(layerName);
            auto coord = parseCoordinate(coordStr);
            startPoints.push_back({layer, coord});
        }
    }
    
    file.close();
    return true;
}

// 点是否在多边形内
bool point_in_polygon(int x, int y, const std::vector<std::pair<int, int>>& vertices) {
    int n = vertices.size();
    bool inside = false;
    
    int p1x = vertices[0].first, p1y = vertices[0].second;
    for (int i = 1; i <= n; i++) {
        int p2x = vertices[i % n].first, p2y = vertices[i % n].second;
        
        if (y > std::min(p1y, p2y)) {
            if (y <= std::max(p1y, p2y)) {
                if (x <= std::max(p1x, p2x)) {
                    double xinters;
                    if (p1y != p2y) {
                        xinters = (double)(y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x;
                    } else {
                        xinters = p1x;
                    }
                    if (p1x == p2x || x <= xinters) {
                        inside = !inside;
                    }
                }
            }
        }
        p1x = p2x;
        p1y = p2y;
    }
    
    return inside;
}

// 查找包含点的多边形
int find_polygon_at(int layer, int x, int y) {
    if (layerIndex.find(layer) == layerIndex.end()) {
        return -1;
    }
    
    Point pt(x, y);
    Box queryBox(pt, pt);
    
    std::vector<RTreeValue> candidates;
    candidates.reserve(16);
    layerIndex[layer].query(bgi::intersects(queryBox), std::back_inserter(candidates));
    
    for (const auto& candidate : candidates) {
        int polyId = candidate.second;
        if (point_in_polygon(x, y, polygons[polyId].vertices)) {
            return polyId;
        }
    }
    
    return -1;
}

// 快速矩形相交判断
inline bool fast_bbox_intersects(const Box& b1, const Box& b2) {
    int minX1 = b1.min_corner().x();
    int maxX1 = b1.max_corner().x();
    int minX2 = b2.min_corner().x();
    int maxX2 = b2.max_corner().x();
    
    if (maxX1 < minX2 || minX1 > maxX2) {
        return false;
    }
    
    int minY1 = b1.min_corner().y();
    int maxY1 = b1.max_corner().y();
    int minY2 = b2.min_corner().y();
    int maxY2 = b2.max_corner().y();
    
    if (maxY1 < minY2 || minY1 > maxY2) {
        return false;
    }
    
    return true;
}

// 判断两个bbox是否接触或重叠
inline bool aabb_touch_or_overlap(const Box& a, const Box& b) {
    return !(b.min_corner().x() > a.max_corner().x() ||
             b.max_corner().x() < a.min_corner().x() ||
             b.min_corner().y() > a.max_corner().y() ||
             b.max_corner().y() < a.min_corner().y());
}

inline Box aabb_intersection(const Box& a, const Box& b) {
    int xmin = std::max(a.min_corner().x(), b.min_corner().x());
    int ymin = std::max(a.min_corner().y(), b.min_corner().y());
    int xmax = std::min(a.max_corner().x(), b.max_corner().x());
    int ymax = std::min(a.max_corner().y(), b.max_corner().y());
    return Box(Point(xmin, ymin), Point(xmax, ymax));
}

enum class ThroughKind { None, VerticalCut, HorizontalCut };

inline ThroughKind through_kind_and_cutC(const PolygonInfo& AA,
                                         const PolygonInfo& PO,
                                         int& c_out) {
    const Box& aaBox = AA.bbox;
    const Box& poBox = PO.bbox;
    if (!aabb_touch_or_overlap(aaBox, poBox)) {
        return ThroughKind::None;
    }

    Box I = aabb_intersection(aaBox, poBox);
    int aa_xmin = aaBox.min_corner().x();
    int aa_xmax = aaBox.max_corner().x();
    int aa_ymin = aaBox.min_corner().y();
    int aa_ymax = aaBox.max_corner().y();
    int ixmin = I.min_corner().x();
    int ixmax = I.max_corner().x();
    int iymin = I.min_corner().y();
    int iymax = I.max_corner().y();

    if (iymin == aa_ymin && iymax == aa_ymax && ixmin < ixmax) {
        long long mid = (static_cast<long long>(ixmin) + static_cast<long long>(ixmax)) >> 1;
        int c = static_cast<int>(mid);
        if (c <= aa_xmin) c = aa_xmin + 1;
        if (c >= aa_xmax) c = aa_xmax - 1;
        if (c > aa_xmin && c < aa_xmax) {
            c_out = c;
            return ThroughKind::VerticalCut;
        }
    }

    if (ixmin == aa_xmin && ixmax == aa_xmax && iymin < iymax) {
        long long mid = (static_cast<long long>(iymin) + static_cast<long long>(iymax)) >> 1;
        int c = static_cast<int>(mid);
        if (c <= aa_ymin) c = aa_ymin + 1;
        if (c >= aa_ymax) c = aa_ymax - 1;
        if (c > aa_ymin && c < aa_ymax) {
            c_out = c;
            return ThroughKind::HorizontalCut;
        }
    }

    return ThroughKind::None;
}

inline bool is_through(const Box& aa, const Box& po) {
    if (!aabb_touch_or_overlap(aa, po)) return false;

    int aa_xmin = aa.min_corner().x();
    int aa_xmax = aa.max_corner().x();
    int aa_ymin = aa.min_corner().y();
    int aa_ymax = aa.max_corner().y();

    int po_xmin = po.min_corner().x();
    int po_xmax = po.max_corner().x();
    int po_ymin = po.min_corner().y();
    int po_ymax = po.max_corner().y();

    bool vertical = (po_xmin <= aa_xmin) && (po_xmax >= aa_xmax) &&
                    !(po_ymax < aa_ymin || po_ymin > aa_ymax);
    bool horizontal = (po_ymin <= aa_ymin) && (po_ymax >= aa_ymax) &&
                      !(po_xmax < aa_xmin || po_xmin > aa_xmax);
    return vertical || horizontal;
}

inline bool gate_stripe_bbox(const PolygonInfo& AA, const PolygonInfo& PO, Box& stripeOut) {
    int cut = 0;
    ThroughKind kind = through_kind_and_cutC(AA, PO, cut);
    if (kind == ThroughKind::None) {
        return false;
    }
    stripeOut = aabb_intersection(AA.bbox, PO.bbox);
    return true;
}

inline bool intervals_overlap_inclusive(int a1, int a2, int b1, int b2) {
    if (a1 > a2) std::swap(a1, a2);
    if (b1 > b2) std::swap(b1, b2);
    return !(a2 < b1 || b2 < a1);
}

inline bool manhattan_intersects(const PolygonInfo& A, const PolygonInfo& B) {
    if (!aabb_touch_or_overlap(A.bbox, B.bbox)) return false;

    for (const auto& he : A.horizEdges) {
        int y = he[0], x1 = he[1], x2 = he[2];
        for (const auto& ve : B.vertEdges) {
            int x = ve[0], y1 = ve[1], y2 = ve[2];
            if (x >= x1 && x <= x2 && y >= y1 && y <= y2) return true;
        }
    }

    for (const auto& he : B.horizEdges) {
        int y = he[0], x1 = he[1], x2 = he[2];
        for (const auto& ve : A.vertEdges) {
            int x = ve[0], y1 = ve[1], y2 = ve[2];
            if (x >= x1 && x <= x2 && y >= y1 && y <= y2) return true;
        }
    }

    for (const auto& ha : A.horizEdges) {
        for (const auto& hb : B.horizEdges) {
            if (ha[0] == hb[0] && intervals_overlap_inclusive(ha[1], ha[2], hb[1], hb[2])) {
                return true;
            }
        }
    }

    for (const auto& va : A.vertEdges) {
        for (const auto& vb : B.vertEdges) {
            if (va[0] == vb[0] && intervals_overlap_inclusive(va[1], va[2], vb[1], vb[2])) {
                return true;
            }
        }
    }

    return false;
}

// 检查两个多边形是否相交（基础版本，无缓存）
inline bool polygons_intersect(int id1, int id2) {
    const auto& p1 = polygons[id1];
    const auto& p2 = polygons[id2];

    if (!fast_bbox_intersects(p1.bbox, p2.bbox)) {
        return false;
    }

    if (g_fastgeo) {
        return manhattan_intersects(p1, p2);
    }
    return bg::intersects(p1.geom, p2.geom);
}

inline bool gate_passable(int aaId, int poId, const std::vector<uint8_t>& isHighPO, Box* stripeOpt = nullptr) {
    if (poId < 0 || poId >= (int)isHighPO.size()) return false;
    if (!isHighPO[poId]) return false;
    if (aaId < 0 || aaId >= (int)polygons.size()) return false;

    const auto& AA = polygons[aaId];
    const auto& PO = polygons[poId];
    if (!polygons_intersect(aaId, poId)) return false;

    Box stripe;
    if (!gate_stripe_bbox(AA, PO, stripe)) return false;
    if (stripeOpt) {
        *stripeOpt = stripe;
    }
    return true;
}

struct SliceTask {
    int aaId;
    std::vector<int> vCuts;
    std::vector<int> hCuts;
};

inline void uniq_sort(std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

std::vector<SliceTask>
collect_slice_tasks(const std::vector<uint8_t>& aaReached,
                    const std::vector<uint8_t>& isHighPO,
                    int aaLayer,
                    int poLayer) {
    std::vector<SliceTask> tasks;
    if (aaLayer < 0 || poLayer < 0) {
        return tasks;
    }

    std::vector<int> poIds;
    if (poLayer < (int)g_layerPolyIds.size()) {
        for (int pid : g_layerPolyIds[poLayer]) {
            if (pid >= 0 && pid < (int)isHighPO.size() && isHighPO[pid]) {
                poIds.push_back(pid);
            }
        }
    } else {
        for (int pid = 0; pid < (int)polygons.size(); ++pid) {
            if (polygons[pid].layer == poLayer && pid < (int)isHighPO.size() && isHighPO[pid]) {
                poIds.push_back(pid);
            }
        }
    }

    if (poIds.empty()) {
        return tasks;
    }

    for (int aaId = 0; aaId < (int)polygons.size(); ++aaId) {
        if (aaId >= (int)aaReached.size() || !aaReached[aaId]) continue;
        const auto& AA = polygons[aaId];
        if (AA.layer != aaLayer) continue;

        SliceTask task;
        task.aaId = aaId;

        for (int poId : poIds) {
            const auto& PO = polygons[poId];
            if (!aabb_touch_or_overlap(AA.bbox, PO.bbox)) continue;
            if (!polygons_intersect(aaId, poId)) continue;

            int c = 0;
            ThroughKind kind = through_kind_and_cutC(AA, PO, c);
            if (kind == ThroughKind::VerticalCut) {
                task.vCuts.push_back(c);
            } else if (kind == ThroughKind::HorizontalCut) {
                task.hCuts.push_back(c);
            }
        }

        auto clip_line_inside = [&](std::vector<int>& cuts, bool vertical) {
            int lo = vertical ? AA.bbox.min_corner().x() : AA.bbox.min_corner().y();
            int hi = vertical ? AA.bbox.max_corner().x() : AA.bbox.max_corner().y();
            std::vector<int> filtered;
            filtered.reserve(cuts.size());
            for (int c : cuts) {
                if (c > lo && c < hi) {
                    filtered.push_back(c);
                }
            }
            cuts.swap(filtered);
            uniq_sort(cuts);
        };

        clip_line_inside(task.vCuts, true);
        clip_line_inside(task.hCuts, false);

        if (!task.vCuts.empty() || !task.hCuts.empty()) {
            tasks.push_back(std::move(task));
        }
    }

    return tasks;
}

inline long long poly_area2(const std::vector<std::pair<int, int>>& V) {
    long long s = 0;
    int n = static_cast<int>(V.size());
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        s += static_cast<long long>(V[i].first) * V[j].second -
             static_cast<long long>(V[j].first) * V[i].second;
    }
    return s;
}

inline bool is_degenerate(const std::vector<std::pair<int, int>>& V) {
    if (static_cast<int>(V.size()) < 3) return true;
    return poly_area2(V) == 0;
}

std::vector<std::pair<int, int>>
clip_vertical(const std::vector<std::pair<int, int>>& V, int c, bool keepLeft) {
    std::vector<std::pair<int, int>> out;
    if (V.empty()) return out;

    auto inside = [&](int x) {
        return keepLeft ? (x <= c) : (x > c);
    };

    int n = static_cast<int>(V.size());
    for (int i = 0; i < n; ++i) {
        auto S = V[i];
        auto E = V[(i + 1) % n];
        bool Sin = inside(S.first);
        bool Ein = inside(E.first);

        if (S.first == E.first && S.first == c) {
            if (keepLeft) {
                if (out.empty() || out.back() != S) out.push_back(S);
                if (out.empty() || out.back() != E) out.push_back(E);
            }
            continue;
        }

        if (Sin && Ein) {
            if (out.empty() || out.back() != E) out.push_back(E);
        } else if (Sin && !Ein) {
            if (S.second == E.second) {
                out.push_back({c, S.second});
            }
        } else if (!Sin && Ein) {
            if (S.second == E.second) {
                out.push_back({c, S.second});
            }
            out.push_back(E);
        }
    }

    if (out.size() >= 2) {
        std::vector<std::pair<int, int>> tmp;
        tmp.reserve(out.size());
        tmp.push_back(out[0]);
        for (size_t i = 1; i < out.size(); ++i) {
            if (out[i] != tmp.back()) {
                tmp.push_back(out[i]);
            }
        }
        out.swap(tmp);
    }

    if (is_degenerate(out)) out.clear();
    return out;
}

std::vector<std::pair<int, int>>
clip_horizontal(const std::vector<std::pair<int, int>>& V, int c, bool keepBottom) {
    std::vector<std::pair<int, int>> out;
    if (V.empty()) return out;

    auto inside = [&](int y) {
        return keepBottom ? (y <= c) : (y > c);
    };

    int n = static_cast<int>(V.size());
    for (int i = 0; i < n; ++i) {
        auto S = V[i];
        auto E = V[(i + 1) % n];
        bool Sin = inside(S.second);
        bool Ein = inside(E.second);

        if (S.second == E.second && S.second == c) {
            if (keepBottom) {
                if (out.empty() || out.back() != S) out.push_back(S);
                if (out.empty() || out.back() != E) out.push_back(E);
            }
            continue;
        }

        if (Sin && Ein) {
            if (out.empty() || out.back() != E) out.push_back(E);
        } else if (Sin && !Ein) {
            if (S.first == E.first) {
                out.push_back({S.first, c});
            }
        } else if (!Sin && Ein) {
            if (S.first == E.first) {
                out.push_back({S.first, c});
            }
            out.push_back(E);
        }
    }

    if (out.size() >= 2) {
        std::vector<std::pair<int, int>> tmp;
        tmp.reserve(out.size());
        tmp.push_back(out[0]);
        for (size_t i = 1; i < out.size(); ++i) {
            if (out[i] != tmp.back()) {
                tmp.push_back(out[i]);
            }
        }
        out.swap(tmp);
    }

    if (is_degenerate(out)) out.clear();
    return out;
}

std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>>
splitByVerticalLine(const std::vector<std::pair<int, int>>& V, int c) {
    auto L = clip_vertical(V, c, true);
    auto R = clip_vertical(V, c, false);
    return {std::move(L), std::move(R)};
}

std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>>
splitByHorizontalLine(const std::vector<std::pair<int, int>>& V, int c) {
    auto B = clip_horizontal(V, c, true);
    auto T = clip_horizontal(V, c, false);
    return {std::move(B), std::move(T)};
}

std::vector<OutPoly>
slice_one_AA(int aaId, const std::vector<int>& vCuts, const std::vector<int>& hCuts) {
    const auto& A = polygons[aaId];
    std::vector<std::vector<std::pair<int, int>>> parts;
    parts.clear();
    parts.emplace_back(A.vertices);

    for (int c : vCuts) {
        std::vector<std::vector<std::pair<int, int>>> next;
        for (auto& P : parts) {
            if (P.empty()) continue;
            auto halves = splitByVerticalLine(P, c);
            if (!halves.first.empty()) next.push_back(std::move(halves.first));
            if (!halves.second.empty()) next.push_back(std::move(halves.second));
        }
        parts.swap(next);
    }

    for (int c : hCuts) {
        std::vector<std::vector<std::pair<int, int>>> next;
        for (auto& P : parts) {
            if (P.empty()) continue;
            auto halves = splitByHorizontalLine(P, c);
            if (!halves.first.empty()) next.push_back(std::move(halves.first));
            if (!halves.second.empty()) next.push_back(std::move(halves.second));
        }
        parts.swap(next);
    }

    std::vector<OutPoly> out;
    out.reserve(parts.size());
    for (auto& P : parts) {
        if (P.size() >= 3 && !is_degenerate(P)) {
            if (!P.empty() && P.front() == P.back()) {
                P.pop_back();
            }
            out.push_back(OutPoly{A.layer, std::move(P)});
        }
    }
    return out;
}

std::vector<OutPoly>
run_slicing_parallel(const std::vector<SliceTask>& tasks) {
    std::vector<std::vector<OutPoly>> perTask(tasks.size());

    #pragma omp parallel for schedule(static) num_threads(g_slice_threads)
    for (int i = 0; i < static_cast<int>(tasks.size()); ++i) {
        const auto& t = tasks[i];
        perTask[i] = slice_one_AA(t.aaId, t.vCuts, t.hCuts);
    }

    std::vector<OutPoly> all;
    size_t total = 0;
    for (const auto& v : perTask) {
        total += v.size();
    }
    all.reserve(total);
    for (auto& v : perTask) {
        all.insert(all.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
    }
    return all;
}

PolygonInfo make_polygoninfo_from_outpoly(const OutPoly& poly) {
    PolygonInfo out;
    out.id = -1;
    out.layer = poly.layer;
    out.originalId = -1;
    out.vertices = poly.verts;
    if (!out.vertices.empty()) {
        int minX = out.vertices[0].first;
        int maxX = out.vertices[0].first;
        int minY = out.vertices[0].second;
        int maxY = out.vertices[0].second;
        out.geom.outer().clear();
        for (const auto& v : out.vertices) {
            minX = std::min(minX, v.first);
            maxX = std::max(maxX, v.first);
            minY = std::min(minY, v.second);
            maxY = std::max(maxY, v.second);
            bg::append(out.geom.outer(), Point(v.first, v.second));
        }
        bg::append(out.geom.outer(), Point(out.vertices.front().first, out.vertices.front().second));
        bg::correct(out.geom);
        out.bbox = Box(Point(minX, minY), Point(maxX, maxY));
    }
    precompute_edges(out);
    return out;
}

inline bool point_in_bbox(const Box& box, int x, int y) {
    return !(x < box.min_corner().x() || x > box.max_corner().x() ||
             y < box.min_corner().y() || y > box.max_corner().y());
}

inline void query_layer_candidates_linear(int layer, const Box& qbox, std::vector<int>& out) {
    if (layer < 0 || layer >= static_cast<int>(g_layerPolyIds.size())) return;

    if (layer < static_cast<int>(g_grid.size())) {
        const auto& G = g_grid[layer];
        if (G.cell > 0 && !G.buckets.empty()) {
            ensure_seen_stamp();
            int mytick = ++g_seen_tick;
            if (g_seen_tick == std::numeric_limits<int>::max()) {
                std::fill(g_seen_stamp.begin(), g_seen_stamp.end(), 0);
                g_seen_tick = 1;
                mytick = ++g_seen_tick;
            }

            long long qxmin = qbox.min_corner().x();
            long long qymin = qbox.min_corner().y();
            long long qxmax = qbox.max_corner().x();
            long long qymax = qbox.max_corner().y();

            long long dx0 = qxmin - G.xmin; if (dx0 < 0) dx0 = 0;
            long long dy0 = qymin - G.ymin; if (dy0 < 0) dy0 = 0;
            long long dx1 = qxmax - G.xmin; if (dx1 < 0) dx1 = 0;
            long long dy1 = qymax - G.ymin; if (dy1 < 0) dy1 = 0;

            int ix0 = static_cast<int>(dx0 / G.cell);
            int iy0 = static_cast<int>(dy0 / G.cell);
            int ix1 = static_cast<int>(dx1 / G.cell);
            int iy1 = static_cast<int>(dy1 / G.cell);

            for (int ix = ix0; ix <= ix1; ++ix) {
                for (int iy = iy0; iy <= iy1; ++iy) {
                    auto it = G.buckets.find(cell_key(ix, iy));
                    if (it == G.buckets.end()) continue;
                    for (int pid : it->second) {
                        if (pid < 0 || pid >= static_cast<int>(polygons.size())) continue;
                        if (g_seen_stamp[pid] == mytick) continue;
                        const auto& box = polygons[pid].bbox;
                        if (box.max_corner().x() < qxmin || qxmax < box.min_corner().x() ||
                            box.max_corner().y() < qymin || qymax < box.min_corner().y()) {
                            continue;
                        }
                        g_seen_stamp[pid] = mytick;
                        out.push_back(pid);
                    }
                }
            }
            return;
        }
    }

    const auto& ids = g_layerPolyIds[layer];
    int qxmin = qbox.min_corner().x();
    int qxmax = qbox.max_corner().x();
    int qymin = qbox.min_corner().y();
    int qymax = qbox.max_corner().y();

    for (int pid : ids) {
        const auto& box = polygons[pid].bbox;
        if (!(box.max_corner().x() < qxmin || qxmax < box.min_corner().x() ||
              box.max_corner().y() < qymin || qymax < box.min_corner().y())) {
            out.push_back(pid);
        }
    }
}

inline void query_layer_candidates(int layer, const Box& qbox, std::vector<int>& out) {
    if (layer >= 0 && layer < static_cast<int>(g_grid.size())) {
        const auto& G = g_grid[layer];
        if (G.cell > 0 && !G.buckets.empty()) {
            query_layer_candidates_linear(layer, qbox, out);
            return;
        }
    }

    auto it = layerIndex.find(layer);
    if (it != layerIndex.end()) {
        std::vector<RTreeValue> candidates;
        candidates.reserve(64);
        it->second.query(bgi::intersects(qbox), std::back_inserter(candidates));
        for (const auto& candidate : candidates) {
            out.push_back(candidate.second);
        }
        return;
    }

    query_layer_candidates_linear(layer, qbox, out);
}

inline void build_per_layer_registry() {
    int maxL = 0;
    for (const auto& poly : polygons) {
        if (poly.layer > maxL) {
            maxL = poly.layer;
        }
    }
    g_layerPolyIds.assign(maxL + 1, {});
    for (int pid = 0; pid < (int)polygons.size(); ++pid) {
        int L = polygons[pid].layer;
        if (L >= 0) {
            if (L >= (int)g_layerPolyIds.size()) {
                g_layerPolyIds.resize(L + 1);
            }
            g_layerPolyIds[L].push_back(pid);
        }
    }
}

inline bool point_in_poly_manhattan(const PolygonInfo& P, int x, int y) {
    const auto& V = P.vertices;
    if (V.empty()) return false;

    bool inside = false;
    for (size_t i = 0, j = V.size() - 1; i < V.size(); j = i++) {
        int xi = V[i].first, yi = V[i].second;
        int xj = V[j].first, yj = V[j].second;

        if ((y == yi && y == yj && ((x - xi) * (x - xj) <= 0)) ||
            (x == xi && x == xj && ((y - yi) * (y - yj) <= 0))) {
            return true;
        }

        bool intersect = ((yi > y) != (yj > y)) &&
                         (x < (int64_t)(xj - xi) * (y - yi) / (int64_t)(yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
}

inline int locate_seed_polygon_linear(int layerHint, int sx, int sy) {
    auto find_in_layer = [&](int layer) -> int {
        if (layer < 0 || layer >= (int)g_layerPolyIds.size()) return -1;
        Box queryBox(Point(sx, sy), Point(sx, sy));
        std::vector<int> cand;
        cand.reserve(16);
        query_layer_candidates(layer, queryBox, cand);

        if (!cand.empty()) {
            for (int pid : cand) {
                const auto& poly = polygons[pid];
                if (!point_in_bbox(poly.bbox, sx, sy)) continue;
                if (point_in_poly_manhattan(poly, sx, sy)) return pid;
            }
        }

        for (int pid : g_layerPolyIds[layer]) {
            const auto& poly = polygons[pid];
            if (!point_in_bbox(poly.bbox, sx, sy)) continue;
            if (point_in_poly_manhattan(poly, sx, sy)) return pid;
        }
        return -1;
    };

    if (layerHint >= 0) {
        int res = find_in_layer(layerHint);
        if (res >= 0) return res;
    }

    for (int layer = 0; layer < (int)g_layerPolyIds.size(); ++layer) {
        if (layer == layerHint) continue;
        int res = find_in_layer(layer);
        if (res >= 0) return res;
    }

    for (size_t idx = 0; idx < polygons.size(); ++idx) {
        const auto& poly = polygons[idx];
        if (!point_in_bbox(poly.bbox, sx, sy)) continue;
        if (point_in_poly_manhattan(poly, sx, sy)) return static_cast<int>(idx);
    }
    return -1;
}

inline void expand_same_layer_on_demand(int layer, int pid, std::vector<int>& out) {
    if (layer < 0) return;
    const auto& P = polygons[pid];
    std::vector<int> cand;
    cand.reserve(64);
    query_layer_candidates(layer, P.bbox, cand);

    for (int qid : cand) {
        if (qid == pid) continue;
        if (polygons_intersect(pid, qid)) {
            out.push_back(qid);
        }
    }
}

inline void expand_adjacent_layers_on_demand(int layer, int pid,
                                             const std::vector<std::vector<int>>& layerNeighbors,
                                             std::vector<std::pair<int,int>>& out) {
    if (layer < 0 || layer >= (int)layerNeighbors.size()) return;
    const auto& P = polygons[pid];
    for (int L2 : layerNeighbors[layer]) {
        if (L2 < 0) continue;
        std::vector<int> cand;
        cand.reserve(32);
        query_layer_candidates(L2, P.bbox, cand);
        for (int qid : cand) {
            if (polygons_intersect(pid, qid)) {
                out.emplace_back(L2, qid);
            }
        }
    }
}

inline void expand_same_layer_s2(int layer, int uid, std::vector<int>& out, int aaLayer) {
    if (layer < 0) return;
    if (aaLayer >= 0 && layer == aaLayer) return;

    const auto& U = polygons[uid];
    std::vector<int> cand;
    cand.reserve(64);
    query_layer_candidates(layer, U.bbox, cand);
    for (int vid : cand) {
        if (vid == uid) continue;
        if (polygons_intersect(uid, vid)) {
            out.push_back(vid);
        }
    }
}

inline void expand_cross_layer_s2(int layer, int uid,
                                  const std::vector<std::vector<int>>& layerNeighbors,
                                  std::vector<std::pair<int,int>>& out,
                                  const std::vector<uint8_t>& isHighPO,
                                  int aaLayer, int poLayer) {
    if (layer < 0 || layer >= (int)layerNeighbors.size()) return;

    const auto& U = polygons[uid];

    if (aaLayer >= 0 && layer == aaLayer && poLayer >= 0) {
        std::vector<int> pocand;
        pocand.reserve(64);
        query_layer_candidates(poLayer, U.bbox, pocand);

        std::vector<Box> stripes;
        stripes.reserve(pocand.size());
        for (int pid : pocand) {
            Box stripe;
            if (gate_passable(uid, pid, isHighPO, &stripe)) {
                stripes.push_back(stripe);
            }
        }
        if (stripes.empty()) {
            return;
        }

        std::vector<int> cand;
        cand.reserve(64);
        for (int L2 : layerNeighbors[layer]) {
            if (L2 < 0) continue;
            cand.clear();
            query_layer_candidates(L2, U.bbox, cand);
            for (int vid : cand) {
                if (vid == uid) continue;
                if (!polygons_intersect(uid, vid)) continue;
                const auto& V = polygons[vid];
                bool ok = false;
                for (const auto& stripe : stripes) {
                    if (aabb_touch_or_overlap(V.bbox, stripe)) {
                        ok = true;
                        break;
                    }
                }
                if (ok) {
                    out.emplace_back(L2, vid);
                }
            }
        }
        return;
    }

    std::vector<int> cand;
    cand.reserve(64);
    for (int L2 : layerNeighbors[layer]) {
        if (L2 < 0) continue;
        cand.clear();
        query_layer_candidates(L2, U.bbox, cand);
        for (int vid : cand) {
            if (vid == uid) continue;
            if (polygons_intersect(uid, vid)) {
                out.emplace_back(L2, vid);
            }
        }
    }
}

struct Seed {
    int x;
    int y;
    int layerHint;
};

struct RingQ {
    std::vector<int> buf;
    size_t head = 0;
    size_t tail = 0;
    size_t mask = 0;
    explicit RingQ(size_t cap = 1u << 14) {
        size_t c = 1;
        while (c < cap) c <<= 1;
        buf.resize(c);
        mask = buf.size() - 1;
    }
    inline bool empty() const { return head == tail; }
    inline size_t size() const { return tail - head; }
    inline void grow() {
        size_t current = size();
        std::vector<int> next(buf.size() << 1);
        for (size_t i = 0; i < current; ++i) {
            next[i] = buf[(head + i) & mask];
        }
        buf.swap(next);
        head = 0;
        tail = current;
        mask = buf.size() - 1;
    }
    inline void push(int v) {
        if (size() == buf.size()) grow();
        buf[tail & mask] = v;
        ++tail;
    }
    inline int pop() {
        int v = buf[head & mask];
        ++head;
        return v;
    }
};

inline void build_layer_neighbor_table(const std::vector<std::vector<int>>& viaRules,
                                       std::vector<std::vector<int>>& layerNeighbors) {
    int maxLayer = 0;
    for (const auto& poly : polygons) {
        maxLayer = std::max(maxLayer, poly.layer);
    }
    for (const auto& rule : viaRules) {
        for (int layer : rule) {
            maxLayer = std::max(maxLayer, layer);
        }
    }
    layerNeighbors.assign(maxLayer + 1, {});

    for (const auto& rule : viaRules) {
        for (size_t i = 0; i + 1 < rule.size(); ++i) {
            int a = rule[i];
            int b = rule[i + 1];
            if (a < 0 || b < 0) continue;
            layerNeighbors[a].push_back(b);
            layerNeighbors[b].push_back(a);
        }
    }

    for (auto& nbrs : layerNeighbors) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
}

bool bfs_from_s1_ondemand(const Seed& s1,
                          int polyLayer,
                          const std::vector<std::vector<int>>& layerNeighbors,
                          std::vector<uint8_t>& isHighPO) {
    PhaseStopwatch _(g_profile ? "[S1]" : nullptr, g_time_s1);
    isHighPO.assign(polygons.size(), 0);
    if (polygons.empty()) return false;

    int startPid = locate_seed_polygon_linear(s1.layerHint, s1.x, s1.y);
    if (startPid < 0) return false;

    std::vector<uint8_t> visitedLocal(polygons.size(), 0);
    RingQ q(polygons.size() + 16);
    q.push(startPid);
    visitedLocal[startPid] = 1;

    auto mark_if_PO = [&](int pid) {
        if (polyLayer != -1 && polygons[pid].layer == polyLayer) {
            isHighPO[pid] = 1;
        }
    };

    mark_if_PO(startPid);

    std::vector<int> sameLayerNeighbors;
    sameLayerNeighbors.reserve(128);
    std::vector<std::pair<int,int>> crossLayerNeighbors;
    crossLayerNeighbors.reserve(64);

    while (!q.empty()) {
        int u = q.pop();
        int Lu = polygons[u].layer;

        sameLayerNeighbors.clear();
        expand_same_layer_on_demand(Lu, u, sameLayerNeighbors);
        for (int v : sameLayerNeighbors) {
            if (!visitedLocal[v]) {
                visitedLocal[v] = 1;
                q.push(v);
                mark_if_PO(v);
            }
        }

        crossLayerNeighbors.clear();
        expand_adjacent_layers_on_demand(Lu, u, layerNeighbors, crossLayerNeighbors);
        for (const auto& item : crossLayerNeighbors) {
            int v = item.second;
            if (!visitedLocal[v]) {
                visitedLocal[v] = 1;
                q.push(v);
                mark_if_PO(v);
            }
        }
    }
    return true;
}

bool bfs_from_s2_with_gate_ondemand(const Seed& s2,
                                    const std::vector<uint8_t>& isHighPO,
                                    const std::vector<std::vector<int>>& layerNeighbors,
                                    int polyLayer, int aaLayer,
                                    std::vector<uint8_t>& aaReached,
                                    std::vector<int>* visitedOut) {
    PhaseStopwatch _(g_profile ? "[S2]" : nullptr, g_time_s2);
    aaReached.assign(polygons.size(), 0);
    if (visitedOut) visitedOut->clear();
    if (polygons.empty()) return false;

    int startPid = locate_seed_polygon_linear(s2.layerHint, s2.x, s2.y);
    if (startPid < 0) return false;

    std::vector<uint8_t> visitedLocal(polygons.size(), 0);
    RingQ q(polygons.size() + 16);
    q.push(startPid);
    visitedLocal[startPid] = 1;

    auto mark_if_AA = [&](int pid) {
        if (aaLayer >= 0 && polygons[pid].layer == aaLayer) {
            aaReached[pid] = 1;
        }
    };

    mark_if_AA(startPid);

    std::vector<int> sameLayerNeighbors;
    sameLayerNeighbors.reserve(128);
    std::vector<std::pair<int,int>> crossLayerNeighbors;
    crossLayerNeighbors.reserve(128);

    while (!q.empty()) {
        int u = q.pop();

        if (visitedOut) {
            visitedOut->push_back(u);
        }

        int Lu = polygons[u].layer;

        sameLayerNeighbors.clear();
        expand_same_layer_s2(Lu, u, sameLayerNeighbors, aaLayer);
        for (int v : sameLayerNeighbors) {
            if (!visitedLocal[v]) {
                visitedLocal[v] = 1;
                q.push(v);
                mark_if_AA(v);
            }
        }

        crossLayerNeighbors.clear();
        expand_cross_layer_s2(Lu, u, layerNeighbors, crossLayerNeighbors,
                              isHighPO, aaLayer, polyLayer);
        for (const auto& item : crossLayerNeighbors) {
            int v = item.second;
            if (!visitedLocal[v]) {
                visitedLocal[v] = 1;
                q.push(v);
                mark_if_AA(v);
            }
        }
    }

    if (aaLayer >= 0) {
        for (size_t i = 0; i < aaReached.size(); ++i) {
            if (aaReached[i] && polygons[i].layer != aaLayer) {
                aaReached[i] = 0;
            }
        }
    } else {
        std::fill(aaReached.begin(), aaReached.end(), 0);
    }

    return true;
}

// 输出结果
bool write_result(const std::string& filename, const std::vector<PolygonInfo>& outputPolygons) {
    PhaseStopwatch _(g_profile ? "[Output]" : nullptr, g_time_output);

    FastWriter writer(std::max<size_t>(outputPolygons.size() * 128ull, 1ull << 20));

    std::map<int, std::vector<const PolygonInfo*>> layerGroups;
    for (const auto& poly : outputPolygons) {
        layerGroups[poly.layer].push_back(&poly);
    }

    for (const auto& entry : layerGroups) {
        int layerId = entry.first;
        const auto& polys = entry.second;
        auto itName = layerIdToName.find(layerId);
        if (itName != layerIdToName.end()) {
            writer.puts(itName->second);
        } else {
            writer.puti(layerId);
        }
        writer.nl();

        for (const auto* poly : polys) {
            const auto& vertices = poly->vertices;
            for (size_t j = 0; j < vertices.size(); ++j) {
                if (j > 0) writer.put(',');
                writer.put('(');
                writer.puti(vertices[j].first);
                writer.put(',');
                writer.puti(vertices[j].second);
                writer.put(')');
            }
            writer.nl();
        }
    }

    if (!writer.write_to_file(filename)) {
        std::cerr << "Error: Cannot open output file: " << filename << "\n";
        return false;
    }
    return true;
}

// ==================== 问题类型识别 ====================
int identify_problem_type(const std::vector<std::pair<int, std::pair<int, int>>>& startPoints,
                          const std::vector<std::vector<int>>& viaRules,
                          int polyLayer, int aaLayer) {
    // 检查是否有Gate规则
    bool hasGate = (polyLayer != -1 && aaLayer != -1);
    if (hasGate) {
        return 3; // 问题3
    }
    
    // 检查起点数量
    if (startPoints.size() != 1) {
        return 3; // 多起点，问题3
    }
    
    // 检查Via规则
    if (viaRules.size() != 1) {
        return 2; // Via规则数量不是1，问题2
    }
    
    // 检查Via规则涉及的层数
    const auto& rule = viaRules[0];
    if (rule.size() != 1) {
        return 2; // Via规则含多层，问题2
    }
    
    // 满足所有条件：单起点、单Via规则、单层、无Gate
    return 1;
}

// ==================== 问题1：单层单种子按需BFS ====================
namespace Problem1 {

// 问题1专用：PairHash和相交缓存（局部作用域）
struct PairHash {
    inline size_t operator()(const std::pair<int, int>& p) const {
        size_t h1 = std::hash<int>{}(p.first);
        size_t h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

// 问题1专用相交判断（带缓存）
inline bool polygons_intersect_cached(int id1, int id2, 
                                      std::unordered_map<std::pair<int,int>, bool, PairHash>& cache) {
    const auto& p1 = polygons[id1];
    const auto& p2 = polygons[id2];
    
    if (!fast_bbox_intersects(p1.bbox, p2.bbox)) {
        return false;
    }
    
    // 规范化键（小ID在前）
    auto key = std::make_pair(std::min(id1, id2), std::max(id1, id2));
    auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }
    
    bool result = bg::intersects(p1.geom, p2.geom);
    
    // 只缓存相交为true的对（节省内存）
    if (result) {
        cache[key] = result;
    }
    
    return result;
}
    
std::vector<int> trace_single_layer_on_demand(int seedPolygon, int layer) {
    std::vector<int> region;
    std::queue<int> queue;
    
    if (visited[seedPolygon]) {
        return region;
    }
    
    std::cerr << "Problem 1: Single-layer on-demand BFS (layer=" << layerIdToName[layer] << ")\n";
    std::cerr << "  Optimizations: Candidate cache + Intersection cache\n";
    
    if (layerIndex.find(layer) == layerIndex.end()) {
        std::cerr << "Error: Layer " << layerIdToName[layer] << " has no R-tree index\n";
        return region;
    }
    
    const RTree& rtree = layerIndex[layer];
    
    // 候选缓存（惰性填充）
    std::unordered_map<int, std::vector<int>> neighbors_cache;
    
    // 相交结果缓存（问题1局部）
    std::unordered_map<std::pair<int,int>, bool, PairHash> intersectionCache;
    
    visited[seedPolygon] = true;
    queue.push(seedPolygon);
    region.push_back(seedPolygon);
    
    long long rtree_queries = 0;
    long long cache_hits = 0;
    long long bbox_checks = 0;
    long long bbox_rejected = 0;
    long long precise_checks = 0;
    long long intersection_cache_hits = 0;
    long long edges_found = 0;
    
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();
        
        const auto& curPoly = polygons[cur];
        std::vector<int> neighbors;
        
        // 检查候选缓存
        auto cacheIt = neighbors_cache.find(cur);
        if (cacheIt != neighbors_cache.end()) {
            neighbors = cacheIt->second;
            cache_hits++;
        } else {
            // 首次查询，使用R-tree获取候选
            std::vector<RTreeValue> candidates;
            rtree.query(bgi::intersects(curPoly.bbox), std::back_inserter(candidates));
            rtree_queries++;
            
            for (const auto& cand : candidates) {
                int candId = cand.second;
                if (candId == cur || polygons[candId].layer != layer) continue;
                
                bbox_checks++;
                if (!fast_bbox_intersects(curPoly.bbox, polygons[candId].bbox)) {
                    bbox_rejected++;
                    continue;
                }
                
                precise_checks++;
                
                // 使用带缓存的相交判断
                int id1 = std::min(cur, candId);
                int id2 = std::max(cur, candId);
                auto key = std::make_pair(id1, id2);
                
                auto it = intersectionCache.find(key);
                if (it != intersectionCache.end()) {
                    intersection_cache_hits++;
                    if (it->second) {
                        neighbors.push_back(candId);
                        edges_found++;
                    }
                    continue;
                }
                
                if (polygons_intersect_cached(cur, candId, intersectionCache)) {
                    neighbors.push_back(candId);
                    edges_found++;
                }
            }
            
            // 缓存此多边形的邻居列表
            neighbors_cache[cur] = neighbors;
        }
        
        // 处理邻居
        for (int neighbor : neighbors) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue.push(neighbor);
                region.push_back(neighbor);
            }
        }
    }
    
    std::cerr << "On-demand BFS statistics:\n";
    std::cerr << "  R-tree queries: " << rtree_queries << "\n";
    std::cerr << "  Candidate cache hits: " << cache_hits << "\n";
    std::cerr << "  AABB checks: " << bbox_checks << " (rejected: " << bbox_rejected << ")\n";
    std::cerr << "  Precise checks: " << precise_checks << "\n";
    std::cerr << "  Intersection cache hits: " << intersection_cache_hits << "\n";
    std::cerr << "  Edges found: " << edges_found << "\n";
    std::cerr << "  Region size: " << region.size() << " polygons\n";
    
    return region;
}

int solve(const std::vector<std::pair<int, std::pair<int, int>>>& startPoints) {
    std::cerr << "\n=== PROBLEM 1: Single-Layer, Single-Seed (On-Demand BFS) ===\n";
    
    // 初始化访问标记
    visited.assign(polygons.size(), false);
    
    std::vector<PolygonInfo> outputPolygons;
    
    // 只有一个起点
    int seedLayer = startPoints[0].first;
    int sx = startPoints[0].second.first;
    int sy = startPoints[0].second.second;
    
    int seedPolygon = find_polygon_at(seedLayer, sx, sy);
    
    if (seedPolygon < 0) {
        std::cerr << "Error: Start point not found\n";
        return 1;
    }
    
    // 按需BFS（不构建全图）
    auto t_bfs = Clock::now();
    std::vector<int> region = trace_single_layer_on_demand(seedPolygon, seedLayer);
    std::cerr << "⏱️  On-demand BFS: " << elapsed_ms(t_bfs) << " ms\n";
    
    // 收集输出
    for (int id : region) {
        if (id >= 0 && id < (int)polygons.size()) {
            outputPolygons.push_back(polygons[id]);
        }
    }
    
    // 输出结果
    if (!write_result(outputFile, outputPolygons)) {
        return 1;
    }

    return 0;
}

} // namespace Problem1

// ==================== 问题2：多层Via，标准图构建+BFS ====================
namespace Problem2 {

// 问题2专用：Frontier-Driven按需扩展（修正版）
std::vector<int> trace_multilayer_frontier(int seedPolygon, const std::vector<std::vector<int>>& viaRules) {
    std::vector<int> region;
    std::queue<int> queue;
    
    if (visited[seedPolygon]) {
        return region;
    }
    
    std::cerr << "Problem 2: Frontier-driven on-demand expansion (corrected)\n";
    
    // 步骤1：预处理Via规则，构建邻层映射（修正：只连接相邻层）
    std::unordered_map<int, std::vector<int>> adjLayers;
    for (const auto& rule : viaRules) {
        for (size_t i = 0; i + 1 < rule.size(); ++i) {  // ✅ 修正：只取相邻元素
            int layer1 = rule[i];
            int layer2 = rule[i + 1];
            
            // 双向添加
            auto& neighbors1 = adjLayers[layer1];
            if (std::find(neighbors1.begin(), neighbors1.end(), layer2) == neighbors1.end()) {
                neighbors1.push_back(layer2);
            }
            
            auto& neighbors2 = adjLayers[layer2];
            if (std::find(neighbors2.begin(), neighbors2.end(), layer1) == neighbors2.end()) {
                neighbors2.push_back(layer1);
            }
        }
    }
    
    std::cerr << "Adjacent layers map (corrected):\n";
    for (const auto& kv : adjLayers) {
        std::cerr << "  " << layerIdToName[kv.first] << " -> ";
        for (int adj : kv.second) {
            std::cerr << layerIdToName[adj] << " ";
        }
        std::cerr << "\n";
    }
    
    // 步骤2：Frontier-driven BFS
    visited[seedPolygon] = true;
    queue.push(seedPolygon);
    region.push_back(seedPolygon);
    
    long long rtree_queries = 0;
    long long same_layer_checks = 0;
    long long cross_layer_checks = 0;
    long long edges_found = 0;
    
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();
        
        const auto& curPoly = polygons[cur];
        int curLayer = curPoly.layer;
        
        // 同层扩展
        if (layerIndex.find(curLayer) != layerIndex.end()) {
            std::vector<RTreeValue> rtreeResults;
            layerIndex[curLayer].query(bgi::intersects(curPoly.bbox), std::back_inserter(rtreeResults));
            rtree_queries++;
            
            for (const auto& res : rtreeResults) {
                int candId = res.second;
                if (candId == cur || polygons[candId].layer != curLayer) continue;
                if (visited[candId]) continue;
                
                same_layer_checks++;
                
                if (!fast_bbox_intersects(curPoly.bbox, polygons[candId].bbox)) {
                    continue;
                }
                
                if (polygons_intersect(cur, candId)) {
                    visited[candId] = true;
                    queue.push(candId);
                    region.push_back(candId);
                    edges_found++;
                }
            }
        }
        
        // 跨层扩展（只查询相邻层）
        if (adjLayers.find(curLayer) != adjLayers.end()) {
            for (int adjLayer : adjLayers[curLayer]) {
                if (layerIndex.find(adjLayer) == layerIndex.end()) continue;
                
                std::vector<RTreeValue> rtreeResults;
                layerIndex[adjLayer].query(bgi::intersects(curPoly.bbox), std::back_inserter(rtreeResults));
                rtree_queries++;
                
                for (const auto& res : rtreeResults) {
                    int candId = res.second;
                    if (polygons[candId].layer != adjLayer) continue;
                    if (visited[candId]) continue;
                    
                    cross_layer_checks++;
                    
                    if (!fast_bbox_intersects(curPoly.bbox, polygons[candId].bbox)) {
                        continue;
                    }
                    
                    if (polygons_intersect(cur, candId)) {
                        visited[candId] = true;
                        queue.push(candId);
                        region.push_back(candId);
                        edges_found++;
                    }
                }
            }
        }
    }
    
    std::cerr << "Frontier-driven BFS statistics:\n";
    std::cerr << "  R-tree queries: " << rtree_queries << "\n";
    std::cerr << "  Same-layer checks: " << same_layer_checks << "\n";
    std::cerr << "  Cross-layer checks: " << cross_layer_checks << "\n";
    std::cerr << "  Edges found: " << edges_found << "\n";
    std::cerr << "  Region size: " << region.size() << " polygons\n";
    
    return region;
}

void build_connectivity_graph(const std::vector<std::vector<int>>& viaRules) {
    adjacency.resize(polygons.size());
    
    int numThreads = get_thread_count();
    
    // 同层连通关系
    std::vector<int> layersToProcess;
    for (const auto& layerEntry : layerPolygons) {
        layersToProcess.push_back(layerEntry.first);
    }
    
    std::cerr << "Building same-layer connections (" << numThreads << " threads)...\n";
    
    std::vector<std::vector<std::vector<int>>> threadLocalAdj(numThreads);
    for (int t = 0; t < numThreads; t++) {
        threadLocalAdj[t].resize(polygons.size());
    }
    
    #pragma omp parallel for schedule(dynamic) if(numThreads > 1)
    for (size_t layerIdx = 0; layerIdx < layersToProcess.size(); layerIdx++) {
        int layer = layersToProcess[layerIdx];
        int tid = omp_get_thread_num();
        auto& localAdj = threadLocalAdj[tid];
        
        const auto& polyIds = layerPolygons[layer];
        
        for (size_t i = 0; i < polyIds.size(); i++) {
            int id1 = polyIds[i];
            const auto& p1 = polygons[id1];
            
            std::vector<RTreeValue> candidates;
            candidates.reserve(64);
            layerIndex[layer].query(bgi::intersects(p1.bbox), std::back_inserter(candidates));
            
            for (const auto& candidate : candidates) {
                int id2 = candidate.second;
                if (id1 >= id2) continue;
                
                const auto& p2 = polygons[id2];
                
                if (!fast_bbox_intersects(p1.bbox, p2.bbox)) {
                    continue;
                }
                
                if (polygons_intersect(id1, id2)) {
                    localAdj[id1].push_back(id2);
                    localAdj[id2].push_back(id1);
                }
            }
        }
    }
    
    std::cerr << "Merging same-layer results...\n";
    for (int t = 0; t < numThreads; t++) {
        for (size_t id = 0; id < threadLocalAdj[t].size(); id++) {
            if (!threadLocalAdj[t][id].empty()) {
                adjacency[id].insert(adjacency[id].end(), 
                                    threadLocalAdj[t][id].begin(), 
                                    threadLocalAdj[t][id].end());
            }
        }
    }
    
    // 跨层Via连通关系
    std::set<std::pair<int,int>> validLayerPairs;
    
    for (const auto& rule : viaRules) {
        for (size_t i = 0; i + 1 < rule.size(); i++) {
            int layer1 = rule[i];
            int layer2 = rule[i + 1];
            
            if (layerPolygons.find(layer1) == layerPolygons.end() ||
                layerPolygons.find(layer2) == layerPolygons.end()) {
                continue;
            }
            
            if (layerPolygons[layer1].empty() || layerPolygons[layer2].empty()) {
                continue;
            }
            
            auto layerPair = std::make_pair(std::min(layer1, layer2), std::max(layer1, layer2));
            validLayerPairs.insert(layerPair);
        }
    }
    
    std::vector<std::pair<int,int>> layerPairVec(validLayerPairs.begin(), validLayerPairs.end());
    std::cerr << "Building cross-layer connections for " << layerPairVec.size() 
              << " layer pairs (" << numThreads << " threads)...\n";
    
    threadLocalAdj.clear();
    threadLocalAdj.resize(numThreads);
    for (int t = 0; t < numThreads; t++) {
        threadLocalAdj[t].resize(polygons.size());
    }
    
    std::vector<int> layerPairCounts(layerPairVec.size(), 0);
    
    #pragma omp parallel for schedule(dynamic) if(numThreads > 1)
    for (size_t pairIdx = 0; pairIdx < layerPairVec.size(); pairIdx++) {
        int tid = omp_get_thread_num();
        auto& localAdj = threadLocalAdj[tid];
        
        int layer1 = layerPairVec[pairIdx].first;
        int layer2 = layerPairVec[pairIdx].second;
        
        int connectionCount = 0;
        
        for (int id1 : layerPolygons[layer1]) {
            const auto& p1 = polygons[id1];
            
            std::vector<RTreeValue> candidates;
            candidates.reserve(32);
            layerIndex[layer2].query(bgi::intersects(p1.bbox), std::back_inserter(candidates));
            
            for (const auto& candidate : candidates) {
                int id2 = candidate.second;
                const auto& p2 = polygons[id2];
                
                if (!fast_bbox_intersects(p1.bbox, p2.bbox)) {
                    continue;
                }
                
                if (polygons_intersect(id1, id2)) {
                    localAdj[id1].push_back(id2);
                    localAdj[id2].push_back(id1);
                    connectionCount++;
                }
            }
        }
        
        layerPairCounts[pairIdx] = connectionCount;
    }
    
    std::cerr << "Merging cross-layer results...\n";
    for (int t = 0; t < numThreads; t++) {
        for (size_t id = 0; id < threadLocalAdj[t].size(); id++) {
            if (!threadLocalAdj[t][id].empty()) {
                adjacency[id].insert(adjacency[id].end(), 
                                    threadLocalAdj[t][id].begin(), 
                                    threadLocalAdj[t][id].end());
            }
        }
    }
    
    std::cerr << "Via connections established:\n";
    for (size_t pairIdx = 0; pairIdx < layerPairVec.size(); pairIdx++) {
        if (layerPairCounts[pairIdx] > 0) {
            std::cerr << "  " << layerIdToName[layerPairVec[pairIdx].first] << "-" 
                      << layerIdToName[layerPairVec[pairIdx].second] << ": " 
                      << layerPairCounts[pairIdx] << " connections\n";
        }
    }
    
    // 去重
    std::cerr << "Deduplicating adjacency lists (" << numThreads << " threads)...\n";
    
    #pragma omp parallel for schedule(dynamic, 1000) if(numThreads > 1)
    for (size_t id = 0; id < adjacency.size(); id++) {
        if (adjacency[id].empty()) continue;
        std::sort(adjacency[id].begin(), adjacency[id].end());
        adjacency[id].erase(std::unique(adjacency[id].begin(), adjacency[id].end()), adjacency[id].end());
    }
    
    std::cerr << "Graph building complete.\n";
}

std::vector<int> trace_connectivity_from_seed(int seedPolygon) {
    std::vector<int> region;
    std::queue<int> queue;
    
    if (visited[seedPolygon]) {
        return region;
    }
    
    visited[seedPolygon] = true;
    queue.push(seedPolygon);
    
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();
        region.push_back(cur);
        
        for (int neighbor : adjacency[cur]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue.push(neighbor);
            }
        }
    }
    
    return region;
}

int solve(const std::vector<std::pair<int, std::pair<int, int>>>& startPoints,
          const std::vector<std::vector<int>>& viaRules) {
    std::cerr << "\n=== PROBLEM 2: Multi-Layer, Frontier-Driven On-Demand (Corrected) ===\n";
    
    // 初始化访问标记
    visited.resize(polygons.size(), false);
    
    std::vector<PolygonInfo> outputPolygons;
    
    // 使用修正后的frontier-driven按需扩展
    for (const auto& sp : startPoints) {
        int seedLayer = sp.first;
        int sx = sp.second.first;
        int sy = sp.second.second;
        
        int seedPolygon = find_polygon_at(seedLayer, sx, sy);
        
        if (seedPolygon < 0) {
            std::cerr << "Warning: Start point not found\n";
            continue;
        }
        
        auto t_bfs = Clock::now();
        std::vector<int> region = trace_multilayer_frontier(seedPolygon, viaRules);
        std::cerr << "⏱️  Frontier-driven BFS: " << elapsed_ms(t_bfs) << " ms\n";
        
        for (int id : region) {
            if (id >= 0 && id < (int)polygons.size()) {
                outputPolygons.push_back(polygons[id]);
            }
        }
    }
    
    // 输出结果
    if (!write_result(outputFile, outputPolygons)) {
        return 1;
    }

    return 0;
}

} // namespace Problem2

// ==================== 问题3：Gate规则处理（完整流程） ====================
namespace Problem3 {

// Gate边映射关系
std::unordered_map<int, std::unordered_set<int>> polyCuts;
std::unordered_map<int, std::vector<int>> origAAFragments;

void build_connectivity_graph(const std::vector<std::vector<int>>& viaRules) {
    // 调用Problem2的图构建函数（复用代码）
    Problem2::build_connectivity_graph(viaRules);
}

void preslice_aa_by_poly(int polyLayer, int aaLayer) {
    if (polyLayer == -1 || aaLayer == -1) return;
    
    int numThreads = get_thread_count();
    std::cerr << "Pre-slicing AA layer by Poly layer (" << numThreads << " threads)...\n";
    
    polyCuts.clear();
    origAAFragments.clear();
    
    if (layerPolygons.find(aaLayer) == layerPolygons.end()) return;
    if (layerPolygons.find(polyLayer) == layerPolygons.end()) return;
    
    const auto& aaPolyIds = layerPolygons[aaLayer];
    
    // 线程本地结果容器
    std::vector<std::vector<int>> threadAaToRemove(numThreads);
    std::vector<std::vector<PolygonInfo>> threadNewAAPieces(numThreads);
    std::vector<std::unordered_map<int, std::unordered_set<int>>> threadPolyCuts(numThreads);
    std::vector<std::unordered_map<int, std::vector<int>>> threadOrigAAFragments(numThreads);
    
    // 并行处理每个AA多边形
    #pragma omp parallel if(numThreads > 1)
    {
        int tid = omp_get_thread_num();
        
        // 线程本地工作缓冲区
        std::vector<Polygon> currentPieces;
        std::vector<Polygon> nextPieces;
        std::vector<Polygon> diffResult;
        std::vector<RTreeValue> candidates;
        std::vector<int> intersectingPolyIds;
        
        currentPieces.reserve(16);
        nextPieces.reserve(16);
        diffResult.reserve(8);
        candidates.reserve(32);
        intersectingPolyIds.reserve(8);
        
        #pragma omp for schedule(dynamic, 10) nowait
        for (size_t i = 0; i < aaPolyIds.size(); i++) {
            int aaId = aaPolyIds[i];
        const auto& aaPoly = polygons[aaId];
        
        candidates.clear();
        intersectingPolyIds.clear();
        
        layerIndex[polyLayer].query(bgi::intersects(aaPoly.bbox), std::back_inserter(candidates));
        
        for (const auto& candidate : candidates) {
            int polyId = candidate.second;
            if (polygons_intersect(aaId, polyId)) {
                intersectingPolyIds.push_back(polyId);
            }
        }
        
        if (intersectingPolyIds.empty()) {
                continue;
        }
        
            // 记录到线程本地容器
            threadAaToRemove[tid].push_back(aaId);
        
        for (int polyId : intersectingPolyIds) {
                threadPolyCuts[tid][polyId].insert(aaId);
        }
        
            currentPieces.clear();
        currentPieces.push_back(aaPoly.geom);
        
        for (int polyId : intersectingPolyIds) {
                nextPieces.clear();
                
            for (const auto& piece : currentPieces) {
                    diffResult.clear();
                bg::difference(piece, polygons[polyId].geom, diffResult);
                
                for (auto& diffPiece : diffResult) {
                    if (!diffPiece.outer().empty()) {
                            nextPieces.push_back(std::move(diffPiece));
                    }
                }
            }
                
                currentPieces.swap(nextPieces);
        }
        
            // 处理生成的碎片（暂时使用临时ID）
        for (auto& piece : currentPieces) {
            if (piece.outer().size() < 3) continue;
            
            bg::correct(piece);
            
            PolygonInfo newPiece;
                newPiece.id = -1;  // 临时ID，稍后重新分配
            newPiece.layer = aaLayer;
            newPiece.originalId = aaId;
            
                newPiece.vertices.reserve(piece.outer().size());
            for (const auto& pt : piece.outer()) {
                newPiece.vertices.push_back({bg::get<0>(pt), bg::get<1>(pt)});
            }
            if (!newPiece.vertices.empty() && 
                newPiece.vertices.front() == newPiece.vertices.back()) {
                newPiece.vertices.pop_back();
            }
            
            newPiece.bbox = bg::return_envelope<Box>(piece);
                newPiece.geom = std::move(piece);
                
                threadNewAAPieces[tid].push_back(std::move(newPiece));
            }
        }
    } // end omp parallel
    
    // 合并线程本地结果
    std::cerr << "Merging results from " << numThreads << " threads...\n";
    
    std::vector<int> aaToRemove;
    std::vector<PolygonInfo> newAAPieces;
    
    for (int t = 0; t < numThreads; t++) {
        aaToRemove.insert(aaToRemove.end(), 
                         threadAaToRemove[t].begin(), 
                         threadAaToRemove[t].end());
        newAAPieces.insert(newAAPieces.end(),
                          threadNewAAPieces[t].begin(),
                          threadNewAAPieces[t].end());
    }
    
    // 重新分配ID并建立origAAFragments映射
    int baseId = polygons.size();
    std::unordered_map<int, std::vector<int>> aaToNewPieces;
    
    for (size_t i = 0; i < newAAPieces.size(); i++) {
        int newId = baseId + i;
        newAAPieces[i].id = newId;
        aaToNewPieces[newAAPieces[i].originalId].push_back(newId);
    }
    
    // 建立origAAFragments
    origAAFragments = aaToNewPieces;
    
    // 合并polyCuts
    for (int t = 0; t < numThreads; t++) {
        for (const auto& kv : threadPolyCuts[t]) {
            int polyId = kv.first;
            for (int aaId : kv.second) {
                polyCuts[polyId].insert(aaId);
            }
        }
    }
    
    std::cerr << "AA polygons to slice: " << aaToRemove.size() << "\n";
    std::cerr << "New AA pieces: " << newAAPieces.size() << "\n";
    
    std::set<int> removeSet(aaToRemove.begin(), aaToRemove.end());
    auto& aaList = layerPolygons[aaLayer];
    aaList.erase(std::remove_if(aaList.begin(), aaList.end(),
        [&removeSet](int id) { return removeSet.count(id) > 0; }), aaList.end());
    
    for (auto& piece : newAAPieces) {
        int newId = piece.id;
        polygons.push_back(piece);
        layerPolygons[aaLayer].push_back(newId);
    }
    
    std::vector<RTreeValue> values;
    values.reserve(layerPolygons[aaLayer].size());
    for (int id : layerPolygons[aaLayer]) {
        values.emplace_back(polygons[id].bbox, id);
    }
    layerIndex[aaLayer] = RTree(values.begin(), values.end());
    
    std::cerr << "Pre-slicing completed\n";
}

void rebuild_layer_connections_incremental(
    int targetLayer, 
    const std::vector<std::vector<int>>& viaRules,
    const std::unordered_set<int>& oldLayerPolygonIds) {
    
    int numThreads = get_thread_count();
    std::cerr << "Incremental rebuild for layer " << layerIdToName[targetLayer] 
              << " (" << numThreads << " threads)...\n";
    
    #pragma omp parallel for schedule(dynamic, 500) if(numThreads > 1)
    for (size_t id = 0; id < adjacency.size(); id++) {
        auto& neighbors = adjacency[id];
        neighbors.erase(
            std::remove_if(neighbors.begin(), neighbors.end(), 
                [&oldLayerPolygonIds](int neighborId) {
                    return oldLayerPolygonIds.count(neighborId) > 0;
                }),
            neighbors.end()
        );
    }
    
    for (int oldId : oldLayerPolygonIds) {
        if (oldId < (int)adjacency.size()) {
            adjacency[oldId].clear();
        }
    }
    
    if (adjacency.size() < polygons.size()) {
        adjacency.resize(polygons.size());
    }
    
    const auto& polyIds = layerPolygons[targetLayer];
    std::vector<std::mutex> mutexes(polygons.size());
    
    #pragma omp parallel for schedule(dynamic, 20) if(numThreads > 1)
    for (size_t i = 0; i < polyIds.size(); i++) {
        int id1 = polyIds[i];
        const auto& p1 = polygons[id1];
        
        std::vector<RTreeValue> candidates;
        candidates.reserve(64);
        layerIndex[targetLayer].query(bgi::intersects(p1.bbox), std::back_inserter(candidates));
        
        for (const auto& candidate : candidates) {
            int id2 = candidate.second;
            if (id1 >= id2) continue;
            
            const auto& p2 = polygons[id2];
            
            if (!fast_bbox_intersects(p1.bbox, p2.bbox)) {
                continue;
            }
            
            if (polygons_intersect(id1, id2)) {
                {
                    std::lock_guard<std::mutex> lock(mutexes[id1]);
                    adjacency[id1].push_back(id2);
                }
                {
                    std::lock_guard<std::mutex> lock(mutexes[id2]);
                    adjacency[id2].push_back(id1);
                }
            }
        }
    }
    
    std::set<std::pair<int,int>> relevantLayerPairs;
    
    for (const auto& rule : viaRules) {
        for (size_t i = 0; i + 1 < rule.size(); i++) {
            int layer1 = rule[i];
            int layer2 = rule[i + 1];
            
            if (layer1 != targetLayer && layer2 != targetLayer) {
                continue;
            }
            
            if (layerPolygons.find(layer1) == layerPolygons.end() ||
                layerPolygons.find(layer2) == layerPolygons.end()) {
                continue;
            }
            
            if (layerPolygons[layer1].empty() || layerPolygons[layer2].empty()) {
                continue;
            }
            
            auto layerPair = std::make_pair(std::min(layer1, layer2), std::max(layer1, layer2));
            relevantLayerPairs.insert(layerPair);
        }
    }
    
    std::vector<std::pair<int,int>> layerPairVec(relevantLayerPairs.begin(), relevantLayerPairs.end());
    std::vector<int> connectionCounts(layerPairVec.size(), 0);
    
    #pragma omp parallel for schedule(dynamic) if(numThreads > 1)
    for (size_t pairIdx = 0; pairIdx < layerPairVec.size(); pairIdx++) {
        int layer1 = layerPairVec[pairIdx].first;
        int layer2 = layerPairVec[pairIdx].second;
        
        int connectionCount = 0;
        
        int queryLayer = (layer1 == targetLayer) ? layer1 : layer2;
        int indexLayer = (layer1 == targetLayer) ? layer2 : layer1;
        
        for (int id1 : layerPolygons[queryLayer]) {
            const auto& p1 = polygons[id1];
            
            std::vector<RTreeValue> candidates;
            candidates.reserve(32);
            layerIndex[indexLayer].query(bgi::intersects(p1.bbox), std::back_inserter(candidates));
            
            for (const auto& candidate : candidates) {
                int id2 = candidate.second;
                const auto& p2 = polygons[id2];
                
                if (!fast_bbox_intersects(p1.bbox, p2.bbox)) {
                    continue;
                }
                
                if (polygons_intersect(id1, id2)) {
                    {
                        std::lock_guard<std::mutex> lock(mutexes[id1]);
                        adjacency[id1].push_back(id2);
                    }
                    {
                        std::lock_guard<std::mutex> lock(mutexes[id2]);
                        adjacency[id2].push_back(id1);
                    }
                    connectionCount++;
                }
            }
        }
        
        connectionCounts[pairIdx] = connectionCount;
    }
    
    std::cerr << "  Deduplicating adjacency lists...\n";
    #pragma omp parallel for schedule(dynamic, 500) if(numThreads > 1)
    for (size_t id = 0; id < adjacency.size(); id++) {
        if (adjacency[id].empty()) continue;
        std::sort(adjacency[id].begin(), adjacency[id].end());
        adjacency[id].erase(std::unique(adjacency[id].begin(), adjacency[id].end()), adjacency[id].end());
    }
    
    std::cerr << "Incremental rebuild complete.\n";
}

void build_gate_edges_fast(int polyLayer, int aaLayer) {
    if (polyLayer == -1 || aaLayer == -1) return;
    
    int numThreads = get_thread_count();
    std::cerr << "Building gate edges (" << numThreads << " threads, using preslice mappings)...\n";
    
    gateAdjacency.assign(polygons.size(), std::vector<GateEdge>());
    
    std::vector<std::pair<int, std::unordered_set<int>>> polyCutsVec(polyCuts.begin(), polyCuts.end());
    std::vector<long long> threadEdgeCounts(numThreads, 0);
    std::vector<std::mutex> gateMutexes(polygons.size());
    
    #pragma omp parallel for schedule(dynamic, 10) if(numThreads > 1)
    for (size_t idx = 0; idx < polyCutsVec.size(); idx++) {
        int tid = omp_get_thread_num();
        int poId = polyCutsVec[idx].first;
        const auto& PO = polygons[poId];
        const auto& cutOrigAAs = polyCutsVec[idx].second;
        
        for (int origAA : cutOrigAAs) {
            auto it = origAAFragments.find(origAA);
            if (it == origAAFragments.end()) continue;
            
            const auto& frags = it->second;
            if (frags.size() < 2) continue;
            
            std::vector<int> touching;
            touching.reserve(frags.size());
            
            for (int fId : frags) {
                const auto& F = polygons[fId];
                if (aabb_touch_or_overlap(F.bbox, PO.bbox)) {
                    touching.push_back(fId);
                }
            }
            
            if (touching.size() < 2) continue;
            
            for (size_t i = 0; i < touching.size(); ++i) {
                for (size_t j = i + 1; j < touching.size(); ++j) {
                    GateEdge edge;
                    edge.to = touching[j];
                    edge.polyId = poId;
                    
                    {
                        std::lock_guard<std::mutex> lock(gateMutexes[touching[i]]);
                        gateAdjacency[touching[i]].push_back(edge);
                    }
                    
                    edge.to = touching[i];
                    {
                        std::lock_guard<std::mutex> lock(gateMutexes[touching[j]]);
                        gateAdjacency[touching[j]].push_back(edge);
                    }
                    
                    threadEdgeCounts[tid] += 2;
                }
            }
        }
    }
    
    long long edgeCount = 0;
    for (int t = 0; t < numThreads; t++) {
        edgeCount += threadEdgeCounts[t];
    }
    
    std::cerr << "Gate edges built (fast): " << edgeCount << "\n";
}

std::vector<int> trace_connectivity_from_seed(int seedPolygon) {
    std::vector<int> region;
    std::queue<int> queue;
    
    if (visited[seedPolygon]) {
        return region;
    }
    
    visited[seedPolygon] = true;
    queue.push(seedPolygon);
    
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();
        region.push_back(cur);
        
        for (int neighbor : adjacency[cur]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue.push(neighbor);
            }
        }
    }
    
    return region;
}

std::vector<int> trace_connectivity_from_seed_with_gate(int seedPolygon,
                                                         const std::vector<uint8_t>& isHighPO) {
    std::vector<int> region;
    std::queue<int> queue;
    
    if (visited[seedPolygon]) {
        return region;
    }
    
    visited[seedPolygon] = true;
    queue.push(seedPolygon);
    
    while (!queue.empty()) {
        int cur = queue.front();
        queue.pop();
        region.push_back(cur);
        
        for (int neighbor : adjacency[cur]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                queue.push(neighbor);
            }
        }
        
        if (cur < (int)gateAdjacency.size()) {
            for (const auto& edge : gateAdjacency[cur]) {
                if (!visited[edge.to] && edge.polyId < (int)isHighPO.size() && isHighPO[edge.polyId]) {
                    visited[edge.to] = true;
                    queue.push(edge.to);
                }
            }
        }
    }
    
    return region;
}

int solve(const std::vector<std::pair<int, std::pair<int, int>>>& startPoints,
          const std::vector<std::vector<int>>& viaRules,
          int polyLayer, int aaLayer) {
    std::cerr << "\n=== PROBLEM 3: Gate Rule Processing (Full Pipeline) ===\n";
    std::cerr << "Gate rule: Poly=" << layerIdToName[polyLayer]
              << " AA=" << layerIdToName[aaLayer] << "\n";

    bool graphBuilt = false;
    std::vector<uint8_t> isHighPO;

    bool twoSeeds = (startPoints.size() == 2);
    bool useOnDemandS2 = twoSeeds && g_s2bfs;

    std::vector<std::vector<int>> layerNeighbors;
    if ((twoSeeds && g_s1bfs) || useOnDemandS2) {
        build_layer_neighbor_table(viaRules, layerNeighbors);
    }

    if (twoSeeds) {
        std::cerr << "Two-seed Gate tracing mode\n";

        int s1Layer = startPoints[0].first;
        int s1x = startPoints[0].second.first;
        int s1y = startPoints[0].second.second;

        if (g_s1bfs) {
            Seed s1Seed{ s1x, s1y, s1Layer };
            if (!bfs_from_s1_ondemand(s1Seed, polyLayer, layerNeighbors, isHighPO)) {
                std::cerr << "Error: Start point 1 not found\n";
                return 1;
            }

            int highPolyCount = 0;
            for (uint8_t flag : isHighPO) {
                if (flag) highPolyCount++;
            }
            std::cerr << "High-level Poly count: " << highPolyCount << "\n";
        } else {
            auto t_build_graph = Clock::now();
            build_connectivity_graph(viaRules);
            std::cerr << "⏱️  Initial graph building: " << elapsed_ms(t_build_graph) << " ms\n";
            graphBuilt = true;

            visited.assign(polygons.size(), false);
            int seed1Polygon = find_polygon_at(s1Layer, s1x, s1y);
            if (seed1Polygon < 0) {
                std::cerr << "Error: Start point 1 not found\n";
                return 1;
            }

            std::vector<int> s1Region = trace_connectivity_from_seed(seed1Polygon);

            isHighPO.assign(polygons.size(), 0);
            int highPolyCount = 0;
            for (int id : s1Region) {
                if (polygons[id].layer == polyLayer) {
                    isHighPO[id] = 1;
                    highPolyCount++;
                }
            }

            std::cerr << "High-level Poly count: " << highPolyCount << "\n";
        }

        visited.assign(polygons.size(), false);
    } else {
        auto t_build_graph = Clock::now();
        build_connectivity_graph(viaRules);
        std::cerr << "⏱️  Initial graph building: " << elapsed_ms(t_build_graph) << " ms\n";
        graphBuilt = true;
    }

    if (!graphBuilt && !useOnDemandS2) {
        auto t_build_graph = Clock::now();
        build_connectivity_graph(viaRules);
        std::cerr << "⏱️  Initial graph building: " << elapsed_ms(t_build_graph) << " ms\n";
        graphBuilt = true;
    }

    std::vector<uint8_t> aaReached;
    std::vector<PolygonInfo> outputPolygons;

    if (twoSeeds && useOnDemandS2) {
        int s2Layer = startPoints[1].first;
        int s2x = startPoints[1].second.first;
        int s2y = startPoints[1].second.second;

        Seed s2Seed{ s2x, s2y, s2Layer };
        std::vector<int> s2Region;
        if (!bfs_from_s2_with_gate_ondemand(s2Seed, isHighPO, layerNeighbors,
                                            polyLayer, aaLayer, aaReached, &s2Region)) {
            std::cerr << "Error: Start point 2 not found\n";
            return 1;
        }
        std::cerr << "S2 region size: " << s2Region.size() << "\n";

        std::unordered_set<int> aaCutSet;
        std::vector<PolygonInfo> slicedPieces;

        if (g_slice_enable && aaLayer >= 0 && polyLayer >= 0) {
            PhaseStopwatch sliceTimer(g_profile ? "[Slice]" : nullptr, g_time_slice);
            auto tasks = collect_slice_tasks(aaReached, isHighPO, aaLayer, polyLayer);
            if (!tasks.empty()) {
                auto slicedAA = run_slicing_parallel(tasks);
                aaCutSet.reserve(tasks.size());
                for (const auto& task : tasks) {
                    aaCutSet.insert(task.aaId);
                }
                slicedPieces.reserve(slicedAA.size());
                for (const auto& piece : slicedAA) {
                    slicedPieces.push_back(make_polygoninfo_from_outpoly(piece));
                }
            }
        }

        for (int id : s2Region) {
            if (id < 0 || id >= (int)polygons.size()) continue;
            const auto& poly = polygons[id];
            if (g_slice_enable && aaLayer >= 0 && polyLayer >= 0 &&
                poly.layer == aaLayer && aaCutSet.find(id) != aaCutSet.end()) {
                continue;
            }
            outputPolygons.push_back(poly);
        }

        if (!slicedPieces.empty()) {
            outputPolygons.insert(outputPolygons.end(),
                                  std::make_move_iterator(slicedPieces.begin()),
                                  std::make_move_iterator(slicedPieces.end()));
        }
    } else {
        if (!useOnDemandS2) {
            std::unordered_set<int> oldAAPolygonIds;
            if (aaLayer != -1) {
                auto itAA = layerPolygons.find(aaLayer);
                if (itAA != layerPolygons.end()) {
                    oldAAPolygonIds.insert(itAA->second.begin(), itAA->second.end());
                }
            }

            preslice_aa_by_poly(polyLayer, aaLayer);

            rebuild_layer_connections_incremental(aaLayer, viaRules, oldAAPolygonIds);

            build_gate_edges_fast(polyLayer, aaLayer);

            visited.resize(polygons.size(), false);
        }

        if (twoSeeds) {
            int s2Layer = startPoints[1].first;
            int s2x = startPoints[1].second.first;
            int s2y = startPoints[1].second.second;

            int seed2Polygon = find_polygon_at(s2Layer, s2x, s2y);
            if (seed2Polygon < 0) {
                std::cerr << "Error: Start point 2 not found\n";
                return 1;
            }

            std::vector<int> s2Region = trace_connectivity_from_seed_with_gate(seed2Polygon, isHighPO);
            std::cerr << "S2 region size: " << s2Region.size() << "\n";

            for (int id : s2Region) {
                if (id >= 0 && id < (int)polygons.size()) {
                    outputPolygons.push_back(polygons[id]);
                }
            }
        } else {
            for (const auto& sp : startPoints) {
                int seedLayer = sp.first;
                int sx = sp.second.first;
                int sy = sp.second.second;

                int seedPolygon = find_polygon_at(seedLayer, sx, sy);

                if (seedPolygon < 0) {
                    std::cerr << "Warning: Start point not found\n";
                    continue;
                }

                std::vector<int> region = trace_connectivity_from_seed(seedPolygon);

                for (int id : region) {
                    if (id >= 0 && id < (int)polygons.size()) {
                        outputPolygons.push_back(polygons[id]);
                    }
                }
            }
        }
    }

    if (!write_result(outputFile, outputPolygons)) {
        return 1;
    }

    return 0;
}

} // namespace Problem3

// ==================== Main函数：问题类型识别与分发 ====================
int main(int argc, char* argv[]) {
    
    // 解析命令行参数
    if (!parseCommandLine(argc, argv)) {
        return 1;
    }

    if (g_profile) {
        g_time_parse = g_time_pre = g_time_s1 = g_time_s2 = g_time_slice = g_time_output = 0.0;
    }

    // 解析规则与版图文件
    std::vector<std::pair<int, std::pair<int, int>>> startPoints;
    std::vector<std::vector<int>> viaRules;
    int polyLayer, aaLayer;

    {
        PhaseStopwatch _(g_profile ? "[Parse]" : nullptr, g_time_parse);
        if (!parse_rule_file(ruleFile, startPoints, viaRules, polyLayer, aaLayer)) {
            return 1;
        }
        if (!parse_layout_file(layoutFile)) {
            return 1;
        }
    }

    std::cerr << "Via rules parsed:\n";
    for (size_t i = 0; i < viaRules.size(); i++) {
        std::cerr << "  Via" << (i+1) << ": ";
        for (size_t j = 0; j < viaRules[i].size(); j++) {
            if (j > 0) std::cerr << " ";
            std::cerr << layerIdToName[viaRules[i][j]];
        }
        std::cerr << "\n";
    }

    {
        PhaseStopwatch _(g_profile ? "[Preprocess]" : nullptr, g_time_pre);
        build_per_layer_registry();
        build_uniform_grid_indices();
    }

    std::cerr << "Polygons by layer after parsing:\n";
    for (const auto& p : layerPolygons) {
        std::cerr << "  " << layerIdToName[p.first] << ": " << p.second.size() << " polygons\n";
    }
    
    // 识别问题类型
    int problemType = identify_problem_type(startPoints, viaRules, polyLayer, aaLayer);
    
    std::cerr << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cerr << "║ PROBLEM TYPE IDENTIFIED: " << problemType << "                              ║\n";
    switch (problemType) {
        case 1:
            std::cerr << "║ Strategy: Single-layer on-demand BFS (skip full graph)  ║\n";
            break;
        case 2:
            std::cerr << "║ Strategy: Frontier-driven on-demand (corrected)         ║\n";
            break;
        case 3:
            std::cerr << "║ Strategy: Full Gate processing (preslice + rebuild)     ║\n";
            break;
    }
    std::cerr << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    // 根据问题类型分发到对应的求解器
    int result = 0;
    
    switch (problemType) {
        case 1:
            result = Problem1::solve(startPoints);
            break;
        
        case 2:
            result = Problem2::solve(startPoints, viaRules);
            break;
        
        case 3:
            result = Problem3::solve(startPoints, viaRules, polyLayer, aaLayer);
            break;
        
        default:
            std::cerr << "Error: Unknown problem type\n";
            return 1;
    }
    
    if (result != 0) {
        std::cerr << "Error: Solver returned error code " << result << "\n";
        return result;
    }

    if (g_profile) {
        double total = g_time_parse + g_time_pre + g_time_s1 + g_time_s2 + g_time_slice + g_time_output;
        std::cerr.setf(std::ios::fixed);
        std::cerr.precision(3);
        std::cerr << "[TIME] parse=" << g_time_parse
                  << "s pre=" << g_time_pre
                  << "s s1=" << g_time_s1
                  << "s s2=" << g_time_s2
                  << "s slice=" << g_time_slice
                  << "s out=" << g_time_output
                  << "s total=" << total << "s\n";
    }

    std::cerr << "Result written to " << outputFile << "\n";

    return 0;
}


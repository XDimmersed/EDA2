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

// 多边形信息结构体
struct PolygonInfo {
    int id;
    int layer;
    std::vector<std::pair<int, int>> vertices;
    Polygon geom;
    Box bbox;
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

// 相交判断缓存（仅用于问题1优化）
// 注意：问题2和3构建全图时会产生大量缓存条目，反而降低性能
// 因此只在问题1的 namespace 内局部启用

// 命令行参数
std::string layoutFile, ruleFile, outputFile;
int threadCount = 1;

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
        }
    }
    
    if (layoutFile.empty() || ruleFile.empty() || outputFile.empty()) {
        std::cerr << "Usage: " << argv[0] 
                  << " -layout <file> -rule <file> -output <file> [-thread <n>]\n";
        return false;
    }
    
    set_omp_threads(threadCount);
    
    std::cerr << "Configuration:\n";
    std::cerr << "  Layout: " << layoutFile << "\n";
    std::cerr << "  Rule: " << ruleFile << "\n";
    std::cerr << "  Output: " << outputFile << "\n";
    std::cerr << "  Threads: " << g_threadCount << "\n";
    
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

// 检查两个多边形是否相交（基础版本，无缓存）
inline bool polygons_intersect(int id1, int id2) {
    const auto& p1 = polygons[id1];
    const auto& p2 = polygons[id2];
    
    if (!fast_bbox_intersects(p1.bbox, p2.bbox)) {
        return false;
    }
    
    return bg::intersects(p1.geom, p2.geom);
}

// 输出结果
bool write_result(const std::string& filename, const std::vector<PolygonInfo>& outputPolygons) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open output file: " << filename << "\n";
        return false;
    }
    
    int numThreads = get_thread_count();
    
    // 按层分组
    std::map<int, std::vector<const PolygonInfo*>> layerGroups;
    for (const auto& poly : outputPolygons) {
        layerGroups[poly.layer].push_back(&poly);
    }
    
    // 将layers转为vector以便并行处理
    std::vector<std::pair<int, std::vector<const PolygonInfo*>>> layerVec(
        layerGroups.begin(), layerGroups.end());
    
    // 并行格式化每个layer的字符串
    std::vector<std::string> layerStrings(layerVec.size());
    
    #pragma omp parallel for schedule(dynamic, 1) if(numThreads > 1)
    for (size_t i = 0; i < layerVec.size(); i++) {
        int layerId = layerVec[i].first;
        const auto& polys = layerVec[i].second;
        
        std::ostringstream oss;
        oss.str().reserve(polys.size() * 100);
        
        oss << layerIdToName[layerId] << "\n";
        
        for (const auto* poly : polys) {
            const auto& vertices = poly->vertices;
            
            for (size_t j = 0; j < vertices.size(); j++) {
                if (j > 0) oss << ",";
                oss << "(" << vertices[j].first << "," << vertices[j].second << ")";
            }
            oss << "\n";
        }
        
        layerStrings[i] = oss.str();
    }
    
    // 串行合并输出（保持layer顺序）
    for (const auto& str : layerStrings) {
        file << str;
    }
    
    file.close();
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
    visited.resize(polygons.size(), false);
    
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
    auto t_output = Clock::now();
    if (!write_result(outputFile, outputPolygons)) {
        return 1;
    }
    std::cerr << "⏱️  Output writing: " << elapsed_ms(t_output) << " ms\n";
    
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
    auto t_output = Clock::now();
    if (!write_result(outputFile, outputPolygons)) {
        return 1;
    }
    std::cerr << "⏱️  Output writing: " << elapsed_ms(t_output) << " ms\n";
    
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
                                                         const std::vector<char>& isHighPoly) {
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
                if (!visited[edge.to] && edge.polyId < (int)isHighPoly.size() && isHighPoly[edge.polyId]) {
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
    
    // 构建初始图
    auto t_build_graph = Clock::now();
    build_connectivity_graph(viaRules);
    std::cerr << "⏱️  Initial graph building: " << elapsed_ms(t_build_graph) << " ms\n";
    
    // 保存预切片前的AA层多边形ID
    std::unordered_set<int> oldAAPolygonIds(
        layerPolygons[aaLayer].begin(),
        layerPolygons[aaLayer].end()
    );
        
        // 预切片
    auto t_preslice = Clock::now();
        preslice_aa_by_poly(polyLayer, aaLayer);
    std::cerr << "⏱️  Pre-slicing AA: " << elapsed_ms(t_preslice) << " ms\n";
    
    // 增量重建
    auto t_rebuild = Clock::now();
    rebuild_layer_connections_incremental(aaLayer, viaRules, oldAAPolygonIds);
    std::cerr << "⏱️  Rebuilding graph (incremental): " << elapsed_ms(t_rebuild) << " ms\n";
        
        // 构建Gate边
    auto t_gate = Clock::now();
    build_gate_edges_fast(polyLayer, aaLayer);
    std::cerr << "⏱️  Building Gate edges: " << elapsed_ms(t_gate) << " ms\n";
    
    // 初始化访问标记
    visited.resize(polygons.size(), false);
    
    std::vector<PolygonInfo> outputPolygons;
    
    // 双种子Gate追踪
    if (startPoints.size() == 2) {
        std::cerr << "Two-seed Gate tracing mode\n";
        
        // 第一个起点：获取高电平Poly集合
        int s1Layer = startPoints[0].first;
        int s1x = startPoints[0].second.first;
        int s1y = startPoints[0].second.second;
        
        int seed1Polygon = find_polygon_at(s1Layer, s1x, s1y);
        if (seed1Polygon < 0) {
            std::cerr << "Error: Start point 1 not found\n";
            return 1;
        }
        
        std::vector<int> s1Region = trace_connectivity_from_seed(seed1Polygon);
        
        std::vector<char> isHighPoly(polygons.size(), 0);
        int highPolyCount = 0;
        
        for (int id : s1Region) {
            if (polygons[id].layer == polyLayer) {
                isHighPoly[id] = 1;
                highPolyCount++;
            }
        }
        
        std::cerr << "High-level Poly count: " << highPolyCount << "\n";
        
        // 重置visited
        visited.assign(polygons.size(), false);
        
        // 第二个起点：Gate-aware BFS
        int s2Layer = startPoints[1].first;
        int s2x = startPoints[1].second.first;
        int s2y = startPoints[1].second.second;
        
        int seed2Polygon = find_polygon_at(s2Layer, s2x, s2y);
        if (seed2Polygon < 0) {
            std::cerr << "Error: Start point 2 not found\n";
            return 1;
        }
        
        auto t_bfs = Clock::now();
        std::vector<int> s2Region = trace_connectivity_from_seed_with_gate(seed2Polygon, isHighPoly);
        std::cerr << "⏱️  Gate-aware BFS: " << elapsed_ms(t_bfs) << " ms\n";
        std::cerr << "S2 region size: " << s2Region.size() << "\n";
        
        for (int id : s2Region) {
            if (id >= 0 && id < (int)polygons.size()) {
                outputPolygons.push_back(polygons[id]);
            }
        }
    } else {
        // 普通模式
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
    
    // 输出结果
    auto t_output = Clock::now();
    if (!write_result(outputFile, outputPolygons)) {
        return 1;
    }
    std::cerr << "⏱️  Output writing: " << elapsed_ms(t_output) << " ms\n";
    
    return 0;
}

} // namespace Problem3

// ==================== Main函数：问题类型识别与分发 ====================
int main(int argc, char* argv[]) {
    auto t_total = Clock::now();
    
    // 解析命令行参数
    if (!parseCommandLine(argc, argv)) {
        return 1;
    }
    
    // 解析规则文件
    auto t_parse_rule = Clock::now();
    std::vector<std::pair<int, std::pair<int, int>>> startPoints;
    std::vector<std::vector<int>> viaRules;
    int polyLayer, aaLayer;
    
    if (!parse_rule_file(ruleFile, startPoints, viaRules, polyLayer, aaLayer)) {
        return 1;
    }
    std::cerr << "⏱️  Rule parsing: " << elapsed_ms(t_parse_rule) << " ms\n";
    
    std::cerr << "Via rules parsed:\n";
    for (size_t i = 0; i < viaRules.size(); i++) {
        std::cerr << "  Via" << (i+1) << ": ";
        for (size_t j = 0; j < viaRules[i].size(); j++) {
            if (j > 0) std::cerr << " ";
            std::cerr << layerIdToName[viaRules[i][j]];
        }
        std::cerr << "\n";
    }
    
    // 解析版图文件
    auto t_parse_layout = Clock::now();
    if (!parse_layout_file(layoutFile)) {
        return 1;
    }
    std::cerr << "⏱️  Layout parsing: " << elapsed_ms(t_parse_layout) << " ms\n";
    
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
    
    std::cerr << "\n⏱️  ========== TOTAL TIME: " << elapsed_ms(t_total) << " ms ==========\n";
    std::cerr << "Result written to " << outputFile << "\n";
    
    return 0;
}


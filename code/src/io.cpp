#include "io.hpp"
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>

static const Poly* resolve_poly(const Layout& L, const GateCtx* gate, u64 gid){
  u32 lid = gid_lid(gid);
  u32 pid = gid_pid(gid);
  if(gate && gate->is_slice_gid(gid)){
    return gate->slice_from_gid(gid);
  }
  if(lid>=L.layers.size()) return nullptr;
  const auto& layer = L.layers[lid];
  if(pid>=layer.polys.size()) return nullptr;
  return &layer.polys[pid];
}

// 多边形几何签名：用于去重
struct PolyKey {
  u64 h1, h2;
  bool operator==(const PolyKey& o) const { return h1 == o.h1 && h2 == o.h2; }
  bool operator<(const PolyKey& o) const {
    if(h1 != o.h1) return h1 < o.h1;
    return h2 < o.h2;
  }
};

struct PolyHasher {
  size_t operator()(const PolyKey& k) const {
    return (size_t)(k.h1 * 1315423911ULL) ^ (size_t)k.h2;
  }
};

// 规范化多边形的哈希（正确姿势）
static PolyKey canonical_hash(const std::vector<Pt>& v) {
  if(v.empty()) return {0, 0};
  
  // 步骤1：找到字典序最小的起点（规范化起点）
  size_t min_idx = 0;
  for(size_t i = 1; i < v.size(); ++i) {
    if(v[i].x < v[min_idx].x || (v[i].x == v[min_idx].x && v[i].y < v[min_idx].y)) {
      min_idx = i;
    }
  }
  
  // 步骤2：从规范起点开始计算哈希
  u64 h1 = 0, h2 = 0;
  size_t n = v.size();
  
  // 包含顶点数信息（降低碰撞）
  h1 = n * 1000000007ULL;
  h2 = n * 1000000009ULL;
  
  // 按规范顺序哈希顶点
  for(size_t i = 0; i < n; ++i) {
    size_t idx = (min_idx + i) % n;
    u64 x = (u64)(v[idx].x + 1000000000);  // 避免负数
    u64 y = (u64)(v[idx].y + 1000000000);
    h1 = h1 * 1315423911ULL + x;
    h2 = h2 * 1000000007ULL + y;
  }
  
  // 步骤3：混入bbox信息（进一步降低碰撞）
  if(n > 0) {
    i32 minx = v[0].x, miny = v[0].y, maxx = v[0].x, maxy = v[0].y;
    for(const auto& p : v) {
      minx = std::min(minx, p.x);
      miny = std::min(miny, p.y);
      maxx = std::max(maxx, p.x);
      maxy = std::max(maxy, p.y);
    }
    u64 bbox_hash = ((u64)(minx + 1000000000) << 32) | (u64)(miny + 1000000000);
    h1 ^= bbox_hash;
    bbox_hash = ((u64)(maxx + 1000000000) << 32) | (u64)(maxy + 1000000000);
    h2 ^= bbox_hash;
  }
  
  return {h1, h2};
}

void write_result(const Layout& L, const std::vector<u64>& reachable, const GateCtx* gate, const std::string& out_path){
  // ============================================================================
  // 去重的正确姿势：两层防线
  // ============================================================================
  
  // 【防线1】运行期：visited (在BFS中已完成，保证同一gid只处理一次)
  //   - 每个gid最多入队一次
  //   - AA懒切已缓存plan+slices，杜绝重复切
  
  // 【防线2】输出期：按层几何去重
  //   - 不同层的几何完全相同是合理的（如不同金属层的重叠走线）
  //   - 必须分层去重，避免跨层误删
  
  // 步骤1：gid去重（理论上BFS已保证unique，这里double-check）
  std::vector<u64> unique_reachable = reachable;
  std::sort(unique_reachable.begin(), unique_reachable.end());
  unique_reachable.erase(std::unique(unique_reachable.begin(), unique_reachable.end()), 
                         unique_reachable.end());
  
  size_t gid_duplicates = reachable.size() - unique_reachable.size();
  std::cerr << "[write_result] Total reachable: " << reachable.size() 
            << ", unique gids: " << unique_reachable.size()
            << ", gid duplicates removed: " << gid_duplicates << "\n";
  
  std::unordered_map<u32, std::vector<u64>> layer_gids;
  layer_gids.reserve(L.layers.size());
  
  // 统计来源（用于调试）
  std::unordered_map<u32, size_t> slice_count, orig_count;
  std::unordered_map<u32, size_t> dropped_by_layer;
  
  // 步骤2：按层几何去重（每层独立维护seen集合）
  std::unordered_map<u32, std::unordered_set<PolyKey, PolyHasher>> seen_by_layer;
  size_t dropped_dup = 0;
  
  for(u64 gid : unique_reachable){
    const Poly* poly = resolve_poly(L, gate, gid);
    if(!poly || poly->v.size()<3) continue;
    
    // 统计来源
    if(gate && gate->is_slice_gid(gid)) {
      slice_count[poly->lid]++;
    } else {
      orig_count[poly->lid]++;
    }
    
    // 按层几何去重
    PolyKey key = canonical_hash(poly->v);
    auto& layer_seen = seen_by_layer[poly->lid];
    
    if(!layer_seen.insert(key).second) {
      dropped_dup++;
      dropped_by_layer[poly->lid]++;
      continue;  // 跳过该层内几何相同的多边形
    }
    
    layer_gids[poly->lid].push_back(gid);
  }
  
  // 输出统计信息（帮助诊断问题）
  std::cerr << "\n=== Output Statistics ===\n";
  for(const auto& layer : L.layers) {
    size_t slices = slice_count[layer.lid];
    size_t origs = orig_count[layer.lid];
    size_t dropped = dropped_by_layer[layer.lid];
    if(slices > 0 || origs > 0 || dropped > 0) {
      size_t total_before_dedup = slices + origs + dropped;
      std::cerr << layer.name << ": slices=" << slices
                << " + original=" << origs
                << " - dropped=" << dropped
                << " = " << (slices + origs)
                << " (去重率=" << (total_before_dedup > 0 ? dropped*100/total_before_dedup : 0) << "%)\n";
    }
  }
  std::cerr << "Total dropped duplicates: " << dropped_dup << "\n";
  std::cerr << "========================\n\n";
  
  for(auto &kv : layer_gids){
    auto &vec = kv.second;
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  }
  
  std::ofstream fout(out_path);
  for(const auto& layer : L.layers){
    auto it = layer_gids.find(layer.lid);
    if(it==layer_gids.end() || it->second.empty()) continue;
    fout << layer.name << "\n";
    for(u64 gid : it->second){
      const Poly* poly = resolve_poly(L, gate, gid);
      if(!poly) continue;
      for(size_t i=0;i<poly->v.size();++i){
        fout << "(" << poly->v[i].x << "," << poly->v[i].y << ")";
        if(i+1<poly->v.size()) fout << ",";
      }
      fout << "\n";
    }
  }
}

// 诊断输入数据的几何重复情况（在任何BFS之前调用）
void analyze_input_duplicates(const Layout& L) {
  std::cerr << "\n" << std::string(80, '=') << "\n";
  std::cerr << "输入数据几何重复分析（BFS前）\n";
  std::cerr << std::string(80, '=') << "\n\n";
  
  size_t total_all = 0;
  size_t unique_all = 0;
  size_t dups_all = 0;
  
  for(const auto& layer : L.layers) {
    // 统计每个几何形状出现的次数
    std::unordered_map<PolyKey, std::vector<u32>, PolyHasher> geo_to_pids;
    
    for(const auto& poly : layer.polys) {
      PolyKey key = canonical_hash(poly.v);
      geo_to_pids[key].push_back(poly.pid);
    }
    
    size_t total = layer.polys.size();
    size_t unique = geo_to_pids.size();
    size_t dups = total - unique;
    
    total_all += total;
    unique_all += unique;
    dups_all += dups;
    
    double dup_rate = total > 0 ? (dups * 100.0 / total) : 0.0;
    
    std::cerr << layer.name << ":\n";
    std::cerr << "  总数: " << total << "\n";
    std::cerr << "  唯一几何: " << unique << "\n";
    std::cerr << "  重复几何: " << dups << " (" << dup_rate << "%)\n";
    
    // 找出重复最多的几何形状（top 3）
    if(dups > 0) {
      std::vector<std::pair<size_t, PolyKey>> counts;
      for(const auto& kv : geo_to_pids) {
        if(kv.second.size() > 1) {
          counts.push_back({kv.second.size(), kv.first});
        }
      }
      
      std::sort(counts.begin(), counts.end(), std::greater<>());
      
      std::cerr << "  重复最多的几何形状:\n";
      for(size_t i = 0; i < std::min(size_t(3), counts.size()); ++i) {
        std::cerr << "    #" << (i+1) << ": 出现" << counts[i].first << "次\n";
      }
    }
    std::cerr << "\n";
  }
  
  std::cerr << std::string(80, '-') << "\n";
  std::cerr << "全局统计:\n";
  std::cerr << "  总多边形数: " << total_all << "\n";
  std::cerr << "  唯一几何数: " << unique_all << "\n";
  std::cerr << "  重复几何数: " << dups_all << " (" 
            << (total_all > 0 ? (dups_all * 100.0 / total_all) : 0.0) << "%)\n";
  std::cerr << std::string(80, '=') << "\n\n";
}

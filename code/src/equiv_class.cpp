#include "equiv_class.hpp"
#include "geom.hpp"
#include <iostream>
#include <algorithm>
#include <iomanip>

// 计算多边形的规范化几何键
PolyKey EquivClassManager::compute_canonical_key(const std::vector<Pt>& v) const {
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

void EquivClassManager::build_classes(const Layout* layout, u32 aa_layer_id, u32 poly_layer_id) {
  L = layout;
  aa_lid = aa_layer_id;
  poly_lid = poly_layer_id;
  
  size_t num_layers = L->layers.size();
  layer_classes.resize(num_layers);
  layer_modes.resize(num_layers, STRICT); // 默认严格去重
  
  std::cerr << "\n================================================================================\n";
  std::cerr << "等价类构建（预扫描）\n";
  std::cerr << "================================================================================\n\n";
  
  size_t total_polys = 0;
  size_t total_unique = 0;
  size_t total_classes = 0;
  
  for(u32 lid = 0; lid < num_layers; ++lid) {
    const auto& layer = L->layers[lid];
    const auto& polys = layer.polys;
    
    if(polys.empty()) continue;
    
    // AA层和Poly层不建等价类
    if(lid == aa_lid || lid == poly_lid) {
      std::cerr << layer.name << " (特殊层，跳过等价类构建)\n";
      layer_modes[lid] = STRICT;
      continue;
    }
    
    // 按几何键聚类
    std::unordered_map<PolyKey, std::vector<u32>, PolyHasher> buckets;
    for(u32 pid = 0; pid < polys.size(); ++pid) {
      PolyKey key = compute_canonical_key(polys[pid].v);
      buckets[key].push_back(pid);
    }
    
    // 统计
    size_t poly_count = polys.size();
    size_t unique_count = buckets.size();
    size_t dup_count = poly_count - unique_count;
    double dup_rate = poly_count > 0 ? (double)dup_count / poly_count : 0.0;
    
    // 决定该层的策略
    DedupMode mode;
    if(dup_rate >= 0.20) {
      mode = NONE;  // 高重复率：不去重，搜索期用等价类
      class_enabled_layers.insert(lid);
    } else if(dup_rate >= 0.05) {
      mode = LIGHT; // 中等重复率：轻度去重
      class_enabled_layers.insert(lid);
    } else {
      mode = STRICT; // 低重复率：严格去重
      // 不启用等价类（每个多边形独立）
    }
    layer_modes[lid] = mode;
    
    // 构建等价类
    u32 cid = 0;
    size_t max_class_size = 0;
    for(auto& [key, members] : buckets) {
      if(mode == STRICT && members.size() == 1) {
        // 严格模式且无重复，不建类
        continue;
      }
      
      u32 rep = members[0];
      ClassInfo cls(lid, std::move(members), rep, key);
      max_class_size = std::max(max_class_size, cls.members.size());
      
      // 建立gid -> cid映射
      if(mode != STRICT) {
        for(u32 pid : cls.members) {
          u64 gid = ((u64)lid << 32) | pid;
          gid_to_cid[gid] = cid;
        }
      }
      
      layer_classes[lid].push_back(std::move(cls));
      cid++;
    }
    
    // 日志
    std::cerr << layer.name << ":\n";
    std::cerr << "  总数: " << poly_count << "\n";
    std::cerr << "  唯一几何: " << unique_count << "\n";
    std::cerr << "  重复几何: " << dup_count << " (" << std::fixed << std::setprecision(1) << (dup_rate * 100) << "%)\n";
    std::cerr << "  等价类数: " << layer_classes[lid].size() << "\n";
    std::cerr << "  最大类: " << max_class_size << " 成员\n";
    std::cerr << "  策略: ";
    switch(mode) {
      case NONE: std::cerr << "NONE (搜索期类化，输出期全复原)\n"; break;
      case LIGHT: std::cerr << "LIGHT (搜索期类化，输出期GID去重)\n"; break;
      case STRICT: std::cerr << "STRICT (几何去重)\n"; break;
    }
    std::cerr << "\n";
    
    total_polys += poly_count;
    total_unique += unique_count;
    total_classes += layer_classes[lid].size();
  }
  
  std::cerr << "--------------------------------------------------------------------------------\n";
  std::cerr << "全局统计:\n";
  std::cerr << "  总多边形数: " << total_polys << "\n";
  std::cerr << "  唯一几何数: " << total_unique << "\n";
  std::cerr << "  重复几何数: " << (total_polys - total_unique) << " (" 
            << std::fixed << std::setprecision(1) << ((total_polys - total_unique) * 100.0 / total_polys) << "%)\n";
  std::cerr << "  等价类总数: " << total_classes << "\n";
  std::cerr << "  启用等价类的层数: " << class_enabled_layers.size() << "\n";
  std::cerr << "  gid->cid 映射数: " << gid_to_cid.size() << "\n";
  std::cerr << "================================================================================\n\n";
}

bool EquivClassManager::try_get_class_id(u64 gid, u32& out_cid) const {
  auto it = gid_to_cid.find(gid);
  if(it != gid_to_cid.end()) {
    out_cid = it->second;
    return true;
  }
  return false;
}

const ClassInfo& EquivClassManager::get_class(u32 lid, u32 cid) const {
  return layer_classes[lid][cid];
}

const std::vector<u32>& EquivClassManager::get_members(u32 lid, u32 cid) const {
  return layer_classes[lid][cid].members;
}

u32 EquivClassManager::get_rep_pid(u32 lid, u32 cid) const {
  return layer_classes[lid][cid].rep;
}

bool EquivClassManager::is_class_layer(u32 lid) const {
  return class_enabled_layers.count(lid) > 0;
}

DedupMode EquivClassManager::get_layer_mode(u32 lid) const {
  if(lid >= layer_modes.size()) return STRICT;
  return layer_modes[lid];
}

void EquivClassManager::print_statistics() const {
  std::cerr << "\n等价类统计:\n";
  for(u32 lid = 0; lid < layer_classes.size(); ++lid) {
    if(layer_classes[lid].empty()) continue;
    std::cerr << "  Layer " << lid << ": " << layer_classes[lid].size() << " 类\n";
  }
}

size_t EquivClassManager::get_class_count(u32 lid) const {
  if(lid >= layer_classes.size()) return 0;
  return layer_classes[lid].size();
}


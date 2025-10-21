#include "fdec.hpp"
#include "geom.hpp"
#include "parser.hpp"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <iostream>

void FDEC::build(const Layout& layout, const Rule& rule, const Config& cfg_){
  L=&layout; R=&rule; cfg=cfg_;
  // per-layer grid
  grids.resize(L->layers.size());
  for(const auto& layer : L->layers){
    grids[layer.lid] = Grid::build(layer.polys, cfg.grid_cell_hint);
  }
  parse_via_to_adjacency();
  gate_ctx.init_empty();
}

void FDEC::parse_via_to_adjacency(){
  adj_layers.assign(L->layers.size(), {});
  std::unordered_map<std::string,u32> n2id = L->name2lid;
  auto add_adj = [&](u32 a,u32 b){
    auto &va = adj_layers[a], &vb = adj_layers[b];
    if(std::find(va.begin(), va.end(), b)==va.end()) va.push_back(b);
    if(std::find(vb.begin(), vb.end(), a)==vb.end()) vb.push_back(a);
  };
  
  std::cerr << "[Via Rules] Total chains: " << R->vias.size() << "\n";
  for(size_t chain_idx = 0; chain_idx < R->vias.size(); ++chain_idx) {
    const auto& chain = R->vias[chain_idx];
    std::cerr << "  Chain " << chain_idx << ": ";
    for(const auto& layer : chain) {
      std::cerr << layer << " ";
    }
    std::cerr << "\n";
  }
  
  int connection_count = 0;
  // Via规则是链式连接：每行中只有相邻的层建立连通关系
  // 多行Via规则的效果是取并集
  for(const auto& chain : R->vias){
    // 只连接链中相邻的两层
    for(size_t i = 1; i < chain.size(); ++i) {
      auto it1 = n2id.find(chain[i-1]);
      auto it2 = n2id.find(chain[i]);
      if(it1 != n2id.end() && it2 != n2id.end()) {
        add_adj(it1->second, it2->second);
        connection_count++;
      }
    }
  }
  std::cerr << "  Total connections added: " << connection_count << "\n";
  
  std::cerr << "[adj_layers built]\n";
  for(u32 lid = 0; lid < adj_layers.size(); ++lid) {
    if(!adj_layers[lid].empty() && lid < L->layers.size()) {
      std::cerr << "  " << L->layers[lid].name << " (lid=" << lid << ") -> ";
      for(u32 nlid : adj_layers[lid]) {
        if(nlid < L->layers.size()) {
          std::cerr << L->layers[nlid].name << "(" << nlid << ") ";
        }
      }
      std::cerr << "\n";
    }
  }
}

u32 FDEC::locate_seed_pid(u32 lid, const Pt& p) const{
  const auto& polys = L->layers[lid].polys;
  std::vector<u32> candidates;
  
  for(const auto& P: polys){
    if(p.x < P.bb.minx || p.x > P.bb.maxx || p.y < P.bb.miny || p.y > P.bb.maxy) continue;
    if(point_in_poly_manhattan(p, P)) {
      candidates.push_back(P.pid);
    }
  }
  
  // 如果没有候选，返回失败
  if(candidates.empty()) return UINT32_MAX;
  
  // 如果只有一个候选，直接返回
  if(candidates.size() == 1) return candidates[0];
  
  // 如果有多个候选，选择面积最大的（通常更稳定）
  u32 best_pid = candidates[0];
  long long best_area = 0;
  for(u32 pid : candidates) {
    const auto& poly = polys[pid];
    long long area = (long long)(poly.bb.maxx - poly.bb.minx) * (poly.bb.maxy - poly.bb.miny);
    if(area > best_area) {
      best_area = area;
      best_pid = pid;
    }
  }
  return best_pid;
}

void FDEC::expand_frontier_no_gate(u64 gid, std::queue<u64>& q,
                                   std::unordered_set<u64>& visited,
                                   bool allow_aa_hook) const{
  const Poly* U_ptr = nullptr;
  if(gate_ctx.is_slice_gid(gid)){
    U_ptr = gate_ctx.slice_from_gid(gid);
    if(!U_ptr) return;
  }else{
    u32 lid=gid_lid(gid), pid=gid_pid(gid);
    if(lid >= L->layers.size() || pid >= L->layers[lid].polys.size()) return;
    U_ptr = &L->layers[lid].polys[pid];
  }
  const Poly& U = *U_ptr;
  u32 lid = U.lid;

  // ---- A) 同层连通（必须考虑：面积重叠/边重合/角点接触均算连通） ----
  // 重要：如果当前gid是AA切片，禁用同层扩展（切片只能通过邻接边连接其他切片）
  bool is_aa_slice = (gate_ctx.is_slice_gid(gid) && lid == gate_ctx.aa_lid);
  
  if(!is_aa_slice) {
    std::vector<u32> cand; cand.reserve(128);
    grids[lid].query(U.bb, cand);
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    for(u32 vpid : cand){
      // Fix: Only skip self when frontier is from original layer pool.
      // Lazy-cut slices may reuse pid=0, so unconditional skip would
      // wrongly suppress checking polys whose real pid happens to be 0.
      if(!gate_ctx.is_slice_gid(gid) && vpid == U.pid) continue;
      const Poly& V = L->layers[lid].polys[vpid];
      if(!bbox_overlap(U.bb, V.bb)) continue;
      if(!poly_intersect_manhattan(U, V)) continue;
      u64 vgid = make_gid(lid, vpid);
      if(!visited.insert(vgid).second) continue;
      q.push(vgid);
    }
  }

  // ---- 跨层（Via 相邻层）按需扩张 ----
  for(u32 nlid : adj_layers[lid]){
    std::vector<u32> cand; cand.reserve(128);
    grids[nlid].query(U.bb, cand);
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    for(u32 npid : cand){
      const Poly& V = L->layers[nlid].polys[npid];
      if(!bbox_overlap(U.bb, V.bb)) continue;
      if(!poly_intersect_manhattan(U, V)) continue;
      
      // C) Gate 流程：如果nlid是AA，根据allow_aa_hook决定是否执行Lazy-Cut
      if(nlid == gate_ctx.aa_lid && gate_ctx.aa_lid != UINT32_MAX) {
        if(!allow_aa_hook) {
          // Phase 1：允许把AA当作普通多边形，用于via扩展
          // 但不触发Gate逻辑（不切割、不入切片队列）
          u64 ngid = make_gid(nlid, npid);
          if(!visited.insert(ngid).second) continue;
          q.push(ngid);
          continue;
        }
        // Phase 2：允许AA懒切
        const_cast<GateCtx&>(gate_ctx).execute_cut(*L, npid);
        gate_ctx.enqueue_entry_slices_from(gid, npid, q, visited, *L);
        continue;
      }
      
      u64 ngid = make_gid(nlid, npid);
      if(!visited.insert(ngid).second) continue;
      q.push(ngid);
    }
  }
}

std::vector<u64> FDEC::trace_no_gate(const std::vector<std::pair<u32,u32>>& seed_lid_pid) const{
  gate_ctx.init_empty();
  std::queue<u64> q;
  std::unordered_set<u64> visited; visited.reserve(1<<24); // 增加容量

  for(auto [lid,pid] : seed_lid_pid){
    if(lid==UINT32_MAX || pid==UINT32_MAX) continue;
    u64 gid = make_gid(lid,pid);
    if(visited.insert(gid).second) q.push(gid);
  }

  std::vector<u64> out; out.reserve(1<<20); // 增加容量
  while(!q.empty()){
    u64 u=q.front(); q.pop();
    out.push_back(u);
    // 第一二问没有Gate规则，allow_aa_hook参数无影响，使用默认值true
    expand_frontier_no_gate(u, q, visited, /*allow_aa_hook=*/true);
  }
  return out;
}

std::vector<u64> FDEC::trace_with_gate() const{
  // ---- C) Gate 两阶段主流程：s1 激活 poly_high → s2 BFS + Lazy-Cut 钩子 ----
  if(R->seeds.size()<2) return {};

  // s1 起点
  auto s1 = R->seeds[0];
  auto it1 = L->name2lid.find(s1.layer);
  if(it1==L->name2lid.end()) return {};
  u32 lid1 = it1->second;
  u32 pid1 = locate_seed_pid(lid1, s1.p);
  if(pid1==UINT32_MAX) return {};

  // Gate 层绑定
  auto itPoly = L->name2lid.find(R->poly_layer_name);
  auto itAA   = L->name2lid.find(R->aa_layer_name);
  if(itPoly==L->name2lid.end() || itAA==L->name2lid.end()) return {};
  gate_ctx.init_empty();
  gate_ctx.poly_lid = itPoly->second;
  gate_ctx.aa_lid   = itAA->second;
  gate_ctx.poly_high.assign(L->layers[gate_ctx.poly_lid].polys.size(), 0);
  gate_ctx.build_candidates(*L, cfg.grid_cell_hint); // 仅建立 poly↔aa 粗筛映射（不切割）
  
  std::cerr << "[Gate Init] poly_layer=" << R->poly_layer_name << " (lid=" << gate_ctx.poly_lid << ")\n";
  std::cerr << "[Gate Init] aa_layer=" << R->aa_layer_name << " (lid=" << gate_ctx.aa_lid << ")\n";
  std::cerr << "[Gate Init] poly_high.size()=" << gate_ctx.poly_high.size() << "\n";

  // 阶段1：s1 BFS（访问到 poly_layer 即置 high）
  {
    std::unordered_set<u64> visited1; visited1.reserve(1<<20);
    std::queue<u64> q1;
    u64 g1 = make_gid(lid1,pid1);
    visited1.insert(g1); q1.push(g1);
    
    int poly_high_count = 0;

    std::unordered_map<u32, int> layer_poly_counts;
    
    while(!q1.empty()){
      u64 u=q1.front(); q1.pop();
      u32 lid=gid_lid(u), pid=gid_pid(u);
      layer_poly_counts[lid]++;
      
      if(lid==gate_ctx.poly_lid) {
        if(poly_high_count < 10) {
          std::cerr << "[Phase1] Marking poly_high[" << pid << "] = 1\n";
        }
        gate_ctx.poly_high[pid]=1;
        poly_high_count++;
      }
      // Phase 1：禁止AA懒切（allow_aa_hook=false）
      expand_frontier_no_gate(u, q1, visited1, /*allow_aa_hook=*/false);
    }
    
    std::cerr << "Phase 1: visited " << visited1.size() << " polygons, marked " << poly_high_count << " poly as high\n";
    std::cerr << "  poly_high.size() = " << gate_ctx.poly_high.size() << "\n";
    std::cerr << "  Visited by layer:\n";
    for(const auto& [lid, count] : layer_poly_counts) {
      if(count > 0 && lid < L->layers.size()) {
        std::cerr << "    " << L->layers[lid].name << " (lid=" << lid << "): " << count << " polys\n";
      }
    }
  }

  // 阶段2：s2 BFS（命中 AA 时执行多次切割，并通过导通边扩展）
  // 关键修复：清空Phase 1期间可能意外生成的切割计划和切片缓存
  gate_ctx.aa_plans.clear();
  gate_ctx.aa_slices.clear();
  std::cerr << "[Phase 2 Start] Cleared AA plans and slices from Phase 1\n";
  
  auto s2 = R->seeds[1];
  auto it2 = L->name2lid.find(s2.layer);
  if(it2==L->name2lid.end()) return {};
  u32 lid2 = it2->second;
  u32 pid2 = locate_seed_pid(lid2, s2.p);
  if(pid2==UINT32_MAX) return {};

  std::unordered_set<u64> visited2; visited2.reserve(1<<24);
  std::queue<u64> q2;
  std::vector<u64> out; out.reserve(1<<20);
  
  int aa_cut_count = 0, total_pieces = 0;
  std::unordered_set<u32> counted_aa_pids;  // 记录已统计过的AA

  u64 g2 = make_gid(lid2,pid2);
  visited2.insert(g2); q2.push(g2);

  while(!q2.empty()){
    u64 u=q2.front(); q2.pop();
    
    // 处理切割后的片段
    if(gate_ctx.is_slice_gid(u)){
      out.push_back(u);
      total_pieces++;
      
      // 统计切割数量：从切片gid获取原始AA的pid
      auto key = gate_ctx.gid_to_key(u);
      if(key.aa_pid != UINT32_MAX && counted_aa_pids.find(key.aa_pid) == counted_aa_pids.end()) {
        auto it_slices = gate_ctx.aa_slices.find(key.aa_pid);
        if(it_slices != gate_ctx.aa_slices.end() && it_slices->second.ready) {
          if(it_slices->second.pieces.size() > 1) {
            aa_cut_count++;
          }
          counted_aa_pids.insert(key.aa_pid);
        }
      }
      
      // 常规扩张（到其他层，AA的内部扩张已在enqueue_entry_slices_from中完成）
      // Phase 2：允许AA懒切（allow_aa_hook=true，这是默认值）
      expand_frontier_no_gate(u, q2, visited2, /*allow_aa_hook=*/true);
      continue;
    }

    // 统计切割数量（对于未切割的AA，仍然使用原来的方法）
    u32 lid = gid_lid(u), pid = gid_pid(u);
    if(lid == gate_ctx.aa_lid && counted_aa_pids.find(pid) == counted_aa_pids.end()) {
      auto it_slices = gate_ctx.aa_slices.find(pid);
      if(it_slices != gate_ctx.aa_slices.end() && it_slices->second.ready) {
        if(it_slices->second.pieces.size() > 1) {
          aa_cut_count++;
        }
        counted_aa_pids.insert(pid);
      }
    }

    // **关键修复**：如果这是一个已经被切割的AA原始多边形，跳过！
    // 切割后的片段已经通过 enqueue_entry_slices_from 加入队列，
    // 不应该再输出原始AA
    if(lid == gate_ctx.aa_lid) {
      auto it_slices = gate_ctx.aa_slices.find(pid);
      if(it_slices != gate_ctx.aa_slices.end() && 
         it_slices->second.ready && 
         it_slices->second.pieces.size() > 1) {
        // 这个AA已经被切割，原始多边形不应输出，只输出片段
        // 但仍需要通过它来扩展到其他层
        expand_frontier_no_gate(u, q2, visited2, /*allow_aa_hook=*/true);
        continue;
      }
    }

    // 其他层的正常处理，或者未被切割的AA
    out.push_back(u);
    expand_frontier_no_gate(u, q2, visited2, /*allow_aa_hook=*/true);
  }
  
  // 统计有多少片段没被访问
  int total_generated_pieces = 0;
  int aa_cut_into_2 = 0, aa_cut_into_3_plus = 0, aa_not_cut = 0;
  int total_vlines = 0, total_hlines = 0;
  for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
    if(slices.ready) {
      int pc = slices.pieces.size();
      total_generated_pieces += pc;
      if(pc == 1) aa_not_cut++;
      else if(pc == 2) aa_cut_into_2++;
      else if(pc > 2) aa_cut_into_3_plus++;
    }
  }
  
  // 统计切割线数量和高/低电平分布
  int high_vlines = 0, low_vlines = 0;
  int high_hlines = 0, low_hlines = 0;
  
  for(const auto& [aa_pid, plan] : gate_ctx.aa_plans) {
    if(plan.computed) {
      for(const auto& vline : plan.vertical_lines) {
        total_vlines++;
        if(vline.is_high) high_vlines++;
        else low_vlines++;
      }
      for(const auto& hline : plan.horizontal_lines) {
        total_hlines++;
        if(hline.is_high) high_hlines++;
        else low_hlines++;
      }
    }
  }
  
  // 统计访问的原始AA和被切割AA
  int aa_total_in_layer = (int)L->layers[gate_ctx.aa_lid].polys.size();
  int aa_visited_original = 0;  // 访问到的原始AA（未被切割）
  int aa_visited_cut = 0;       // 访问到的被切割AA（至少一个片段被访问）
  
  for(u64 gid : out) {
    u32 lid = gid_lid(gid);
    if(lid == gate_ctx.aa_lid) {
      aa_visited_original++;
    }
  }
  
  for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
    if(!slices.ready || slices.pieces.size() <= 1) continue;
    
    // 检查是否有任何片段被访问
    bool any_visited = false;
    for(size_t i = 0; i < slices.pieces.size(); ++i) {
      u64 slice_gid = gate_ctx.make_slice_gid(aa_pid, (u32)i);
      if(visited2.count(slice_gid)) {
        any_visited = true;
        break;
      }
    }
    if(any_visited) {
      aa_visited_cut++;
    }
  }
  
  std::cerr << "Phase 2: cut " << aa_cut_count << " AA polygons\n";
  std::cerr << "  Not cut (1 piece): " << aa_not_cut << "\n";
  std::cerr << "  Cut into 2 pieces: " << aa_cut_into_2 << "\n";
  std::cerr << "  Cut into 3+ pieces: " << aa_cut_into_3_plus << "\n";
  std::cerr << "  Total vertical cut lines: " << total_vlines << " (high=" << high_vlines 
            << ", low=" << low_vlines << ")\n";
  std::cerr << "  Total horizontal cut lines: " << total_hlines << " (high=" << high_hlines 
            << ", low=" << low_hlines << ")\n";
  std::cerr << "  Total generated pieces: " << total_generated_pieces << "\n";
  std::cerr << "  Total visited pieces: " << total_pieces << "\n";
  std::cerr << "  Unvisited pieces: " << (total_generated_pieces - total_pieces) << "\n";
  // 分析未访问片段
  int unvisited_isolated = 0;  // 孤立片段（没有高电平邻居）
  int unvisited_blocked = 0;   // 被阻断片段（有邻居但都未访问）
  
  for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
    if(!slices.ready || slices.pieces.size() <= 1) continue;
    
    for(size_t i = 0; i < slices.pieces.size(); ++i) {
      u64 slice_gid = gate_ctx.make_slice_gid(aa_pid, i);
      if(visited2.count(slice_gid)) continue;  // 已访问
      
      // 未访问：检查原因
      if(slices.adjacency[i].empty()) {
        unvisited_isolated++;  // 没有高电平邻居
      } else {
        unvisited_blocked++;   // 有邻居但都未访问
      }
    }
  }
  
  std::cerr << "AA Layer Statistics:\n";
  std::cerr << "  Total AA in layer: " << aa_total_in_layer << "\n";
  std::cerr << "  AA visited (original, not cut): " << aa_visited_original << "\n";
  std::cerr << "  AA visited (cut, >=1 slice): " << aa_visited_cut << "\n";
  std::cerr << "  Expected L1 output: " << (aa_visited_original + total_pieces) << "\n";
  std::cerr << "Unvisited Slice Analysis:\n";
  std::cerr << "  Isolated (no high neighbors): " << unvisited_isolated << "\n";
  std::cerr << "  Blocked (has neighbors but unreached): " << unvisited_blocked << "\n";
  
  // 分析孤立片段的空间分布
  std::vector<int> isolated_x_coords;
  std::vector<int> visited_x_coords;
  
  for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
    if(!slices.ready || slices.pieces.size() <= 1) continue;
    
    for(size_t i = 0; i < slices.pieces.size(); ++i) {
      u64 slice_gid = gate_ctx.make_slice_gid(aa_pid, i);
      const auto& piece = slices.pieces[i];
      int center_x = (piece.bb.minx + piece.bb.maxx) / 2;
      
      if(visited2.count(slice_gid)) {
        visited_x_coords.push_back(center_x);
      } else if(slices.adjacency[i].empty()) {
        // 孤立片段
        isolated_x_coords.push_back(center_x);
      }
    }
  }
  
  // 统计x坐标分布
  auto analyze_x_distribution = [](const std::vector<int>& coords, const std::string& name) {
    if(coords.empty()) {
      std::cerr << name << ": 无数据\n";
      return;
    }
    
    int minx = *std::min_element(coords.begin(), coords.end());
    int maxx = *std::max_element(coords.begin(), coords.end());
    
    // 统计不同区域的数量
    int count_left = 0;   // x < 7010 (低电平PO独占区域)
    int count_middle = 0; // 7010 <= x < 333855 (重叠区域)
    int count_right = 0;  // x >= 333855 (低电平PO独占区域)
    
    for(int x : coords) {
      if(x < 7010) count_left++;
      else if(x < 333855) count_middle++;
      else count_right++;
    }
    
    std::cerr << name << ":\n";
    std::cerr << "  总数: " << coords.size() << "\n";
    std::cerr << "  x范围: [" << minx << ", " << maxx << "]\n";
    std::cerr << "  x<7010 (低电平区): " << count_left << " (" 
              << (100.0*count_left/coords.size()) << "%)\n";
    std::cerr << "  7010≤x<333855 (重叠区): " << count_middle << " (" 
              << (100.0*count_middle/coords.size()) << "%)\n";
    std::cerr << "  x≥333855 (低电平区): " << count_right << " (" 
              << (100.0*count_right/coords.size()) << "%)\n";
  };
  
  std::cerr << "\n=== 片段空间分布分析 ===\n";
  analyze_x_distribution(isolated_x_coords, "孤立片段");
  analyze_x_distribution(visited_x_coords, "已访问片段");
  
  // 抽样检查孤立片段的详细信息
  std::cerr << "\n=== 抽样检查孤立片段（前5个）===\n";
  int sample_count = 0;
  for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
    if(!slices.ready || slices.pieces.size() <= 1) continue;
    if(sample_count >= 5) break;
    
    for(size_t i = 0; i < slices.pieces.size(); ++i) {
      u64 slice_gid = gate_ctx.make_slice_gid(aa_pid, i);
      if(visited2.count(slice_gid)) continue;  // 已访问
      if(!slices.adjacency[i].empty()) continue;  // 不是孤立
      
      // 找到一个孤立片段
      const auto& piece = slices.pieces[i];
      const auto& aa = L->layers[gate_ctx.aa_lid].polys[aa_pid];
      
      std::cerr << "\n孤立片段 #" << sample_count << ":\n";
      std::cerr << "  AA pid=" << aa_pid << ", slice_idx=" << i << "\n";
      std::cerr << "  AA bbox: [" << aa.bb.minx << "," << aa.bb.miny << "] 到 [" 
                << aa.bb.maxx << "," << aa.bb.maxy << "]\n";
      std::cerr << "  片段bbox: [" << piece.bb.minx << "," << piece.bb.miny << "] 到 [" 
                << piece.bb.maxx << "," << piece.bb.maxy << "]\n";
      std::cerr << "  片段中心x: " << (piece.bb.minx + piece.bb.maxx)/2 << "\n";
      
      // 检查这个AA的切割计划
      auto it_plan = gate_ctx.aa_plans.find(aa_pid);
      if(it_plan != gate_ctx.aa_plans.end() && it_plan->second.computed) {
        const auto& plan = it_plan->second;
        int high_v = 0, low_v = 0, high_h = 0, low_h = 0;
        for(const auto& vline : plan.vertical_lines) {
          if(vline.is_high) high_v++; else low_v++;
        }
        for(const auto& hline : plan.horizontal_lines) {
          if(hline.is_high) high_h++; else low_h++;
        }
        std::cerr << "  切割线: vert(high=" << high_v << ", low=" << low_v 
                  << "), horiz(high=" << high_h << ", low=" << low_h << ")\n";
        
        // 打印高电平切线的详细信息
        if(high_v + high_h > 0) {
          std::cerr << "  ⚠ 有高电平切线但片段仍孤立！检查邻接逻辑\n";
          for(const auto& vline : plan.vertical_lines) {
            if(vline.is_high) {
              std::cerr << "    高电平竖线: x=" << vline.coord/2 << ", y∈[" 
                        << vline.span_min2/2 << "," << vline.span_max2/2 << "]\n";
            }
          }
          for(const auto& hline : plan.horizontal_lines) {
            if(hline.is_high) {
              std::cerr << "    高电平横线: y=" << hline.coord/2 << ", x∈[" 
                        << hline.span_min2/2 << "," << hline.span_max2/2 << "]\n";
            }
          }
        } else {
          std::cerr << "  ✓ 没有高电平切线，孤立是合理的\n";
        }
      } else {
        std::cerr << "  ⚠ 没有切割计划（不应该发生）\n";
      }
      
      sample_count++;
      if(sample_count >= 5) break;
    }
  }
  
  std::cerr << "\n分析了 " << sample_count << " 个孤立片段样本\n";
  
  // 分析未切割的AA
  std::cerr << "\n=== 分析未切割的AA ===\n";
  int total_aa = L->layers[gate_ctx.aa_lid].polys.size();
  int cut_aa = 0;
  int uncut_aa_with_po = 0;  // 有PO相交但未切割
  int uncut_aa_no_po = 0;    // 没有PO相交
  
  std::vector<u32> sample_uncut_with_po;  // 保存一些样本用于详细分析
  
  for(u32 aa_pid = 0; aa_pid < total_aa; aa_pid++) {
    auto it_slices = gate_ctx.aa_slices.find(aa_pid);
    bool is_cut = (it_slices != gate_ctx.aa_slices.end() && 
                   it_slices->second.ready && 
                   it_slices->second.pieces.size() > 1);
    
    if(is_cut) {
      cut_aa++;
    } else {
      // 未切割，检查是否有PO相交
      const auto& aa_to_poly = gate_ctx.aa_to_poly_cands;
      if(aa_pid < aa_to_poly.size() && !aa_to_poly[aa_pid].empty()) {
        uncut_aa_with_po++;
        if(sample_uncut_with_po.size() < 5) {
          sample_uncut_with_po.push_back(aa_pid);
        }
      } else {
        uncut_aa_no_po++;
      }
    }
  }
  
  std::cerr << "总AA数: " << total_aa << "\n";
  std::cerr << "已切割: " << cut_aa << " (" << (100.0*cut_aa/total_aa) << "%)\n";
  std::cerr << "未切割(有PO相交): " << uncut_aa_with_po << " (" << (100.0*uncut_aa_with_po/total_aa) << "%) ⚠\n";
  std::cerr << "未切割(无PO相交): " << uncut_aa_no_po << " (" << (100.0*uncut_aa_no_po/total_aa) << "%)\n";
  
  // 抽样检查未切割但有PO相交的AA
  std::cerr << "\n=== 抽样检查未切割但有PO的AA（前5个）===\n";
  for(size_t idx = 0; idx < sample_uncut_with_po.size(); idx++) {
    u32 aa_pid = sample_uncut_with_po[idx];
    const auto& aa = L->layers[gate_ctx.aa_lid].polys[aa_pid];
    const auto& po_cands = gate_ctx.aa_to_poly_cands[aa_pid];
    
    std::cerr << "\n未切割AA #" << idx << ":\n";
    std::cerr << "  AA pid=" << aa_pid << "\n";
    std::cerr << "  AA bbox: [" << aa.bb.minx << "," << aa.bb.miny << "] 到 [" 
              << aa.bb.maxx << "," << aa.bb.maxy << "]\n";
    std::cerr << "  相交PO数量: " << po_cands.size() << "\n";
    
    // 检查这些PO的电平
    int high_po = 0, low_po = 0;
    for(u32 po_pid : po_cands) {
      if(po_pid < gate_ctx.poly_high.size() && gate_ctx.poly_high[po_pid]) {
        high_po++;
      } else {
        low_po++;
      }
    }
    std::cerr << "  高电平PO: " << high_po << ", 低电平PO: " << low_po << "\n";
    
    // 检查贯穿判定
    if(high_po > 0 || low_po > 0) {
      std::cerr << "  检查贯穿判定:\n";
      int sample_count_po = 0;
      for(u32 po_pid : po_cands) {
        if(sample_count_po >= 3) break;  // 只看前3个PO
        
        const auto& po = L->layers[gate_ctx.poly_lid].polys[po_pid];
        bool is_high = (po_pid < gate_ctx.poly_high.size() && gate_ctx.poly_high[po_pid]);
        
        // 计算相交区域
        BBox intersection = {
          std::max(po.bb.minx, aa.bb.minx),
          std::max(po.bb.miny, aa.bb.miny),
          std::min(po.bb.maxx, aa.bb.maxx),
          std::min(po.bb.maxy, aa.bb.maxy)
        };
        
        int gate_width = po.bb.maxx - po.bb.minx;
        int gate_height = po.bb.maxy - po.bb.miny;
        
        std::cerr << "    PO#" << po_pid << " (" << (is_high ? "高" : "低") << "电平):\n";
        std::cerr << "      Gate尺寸: " << gate_width << " x " << gate_height << "\n";
        std::cerr << "      交集: [" << intersection.minx << "," << intersection.miny 
                  << "] 到 [" << intersection.maxx << "," << intersection.maxy << "]\n";
        std::cerr << "      交集尺寸: " << (intersection.maxx-intersection.minx) 
                  << " x " << (intersection.maxy-intersection.miny) << "\n";
        
        // 简单判定
        if(gate_height > gate_width) {
          std::cerr << "      判定方向: 竖切 (gate_height > gate_width)\n";
        } else {
          std::cerr << "      判定方向: 横切 (gate_height <= gate_width)\n";
        }
        
        sample_count_po++;
      }
      
      std::cerr << "  ⚠ 有" << (high_po + low_po) << "个PO相交但AA未被切割！\n";
    }
  }
  
  std::cerr << "\n分析了 " << sample_uncut_with_po.size() << " 个未切割AA样本\n";
  
  // 统计孤立片段的切线情况（按片段级别，不是AA级别）
  std::cerr << "\n=== 统计孤立片段的切线配置 ===\n";
  std::cerr << "注意：统计的是\"该片段周围\"是否有高电平切线，而非整个AA\n";
  int isolated_no_high_nearby = 0;      // 周围没有高电平切线（合理孤立）
  int isolated_has_high_nearby = 0;     // 周围有高电平切线但仍孤立（可能是Bug）
  
  for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
    if(!slices.ready || slices.pieces.size() <= 1) continue;
    
    // 检查这个AA的切割计划
    auto it_plan = gate_ctx.aa_plans.find(aa_pid);
    if(it_plan == gate_ctx.aa_plans.end() || !it_plan->second.computed) continue;
    
    const auto& plan = it_plan->second;
    
    // 检查每个片段
    for(size_t i = 0; i < slices.pieces.size(); ++i) {
      u64 slice_gid = gate_ctx.make_slice_gid(aa_pid, i);
      if(visited2.count(slice_gid)) continue;  // 已访问
      if(!slices.adjacency[i].empty()) continue;  // 不是孤立
      
      // 这是一个孤立片段，检查它周围是否有高电平切线
      const auto& piece_bb = slices.pieces[i].bb;
      bool has_high_line_nearby = false;
      const int proximity = 20;  // 20单位内算"周围"
      
      // 检查竖线
      for(const auto& vline : plan.vertical_lines) {
        if(!vline.is_high) continue;
        int line_x = vline.coord / 2;
        int line_y_min = vline.span_min2 / 2;
        int line_y_max = vline.span_max2 / 2;
        
        // 切线是否在片段附近？
        bool x_near = (line_x >= piece_bb.minx - proximity && line_x <= piece_bb.maxx + proximity);
        bool y_overlap = !(line_y_max < piece_bb.miny || line_y_min > piece_bb.maxy);
        
        if(x_near && y_overlap) {
          has_high_line_nearby = true;
          break;
        }
      }
      
      // 检查横线
      if(!has_high_line_nearby) {
        for(const auto& hline : plan.horizontal_lines) {
          if(!hline.is_high) continue;
          int line_y = hline.coord / 2;
          int line_x_min = hline.span_min2 / 2;
          int line_x_max = hline.span_max2 / 2;
          
          // 切线是否在片段附近？
          bool y_near = (line_y >= piece_bb.miny - proximity && line_y <= piece_bb.maxy + proximity);
          bool x_overlap = !(line_x_max < piece_bb.minx || line_x_min > piece_bb.maxx);
          
          if(y_near && x_overlap) {
            has_high_line_nearby = true;
            break;
          }
        }
      }
      
      if(has_high_line_nearby) {
        isolated_has_high_nearby++;
      } else {
        isolated_no_high_nearby++;
      }
    }
  }
  
  std::cerr << "孤立片段总数: " << (unvisited_isolated + unvisited_blocked) << "\n";
  std::cerr << "  周围有高电平切线: " << isolated_has_high_nearby 
            << " (" << (100.0*isolated_has_high_nearby/(unvisited_isolated+unvisited_blocked)) << "%) ⚠\n";
  std::cerr << "  周围无高电平切线: " << isolated_no_high_nearby 
            << " (" << (100.0*isolated_no_high_nearby/(unvisited_isolated+unvisited_blocked)) << "%) ✓\n";
  
  if(isolated_has_high_nearby > 0) {
    std::cerr << "\n⚠ 发现 " << isolated_has_high_nearby << " 个片段周围有高电平切线但仍孤立！\n";
    std::cerr << "   这可能是邻接判定的Bug，或者切割方向错误导致。\n";
    
    // 详细检查前2个这样的孤立片段
    std::cerr << "\n=== 详细检查错误孤立片段（前2个）===\n";
    int debug_count = 0;
    for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
      if(debug_count >= 2) break;
      if(!slices.ready || slices.pieces.size() <= 1) continue;
      
      auto it_plan = gate_ctx.aa_plans.find(aa_pid);
      if(it_plan == gate_ctx.aa_plans.end() || !it_plan->second.computed) continue;
      
      const auto& plan = it_plan->second;
      bool has_high_line = false;
      for(const auto& vline : plan.vertical_lines) {
        if(vline.is_high) { has_high_line = true; break; }
      }
      for(const auto& hline : plan.horizontal_lines) {
        if(hline.is_high) { has_high_line = true; break; }
      }
      if(!has_high_line) continue;
      
      // 找到一个孤立片段
      for(size_t i = 0; i < slices.pieces.size(); ++i) {
        u64 slice_gid = gate_ctx.make_slice_gid(aa_pid, i);
        if(visited2.count(slice_gid)) continue;
        if(!slices.adjacency[i].empty()) continue;
        
        // 这是一个"有高电平切线但孤立"的片段
        const auto& piece = slices.pieces[i];
        const auto& aa = L->layers[gate_ctx.aa_lid].polys[aa_pid];
        
        std::cerr << "\n孤立片段案例 #" << debug_count << ":\n";
        std::cerr << "  AA pid=" << aa_pid << ", 片段idx=" << i << "/" << slices.pieces.size() << "\n";
        std::cerr << "  AA bbox: [" << aa.bb.minx << "," << aa.bb.miny << "] 到 [" 
                  << aa.bb.maxx << "," << aa.bb.maxy << "]\n";
        std::cerr << "  孤立片段bbox: [" << piece.bb.minx << "," << piece.bb.miny << "] 到 [" 
                  << piece.bb.maxx << "," << piece.bb.maxy << "]\n";
        
        // 打印所有片段和切线
        std::cerr << "  所有片段:\n";
        for(size_t j = 0; j < slices.pieces.size(); j++) {
          const auto& p = slices.pieces[j];
          std::cerr << "    #" << j << ": [" << p.bb.minx << "," << p.bb.miny << "] 到 [" 
                    << p.bb.maxx << "," << p.bb.maxy << "], 邻居=" << slices.adjacency[j].size();
          if(j == i) std::cerr << " ← 孤立片段";
          std::cerr << "\n";
        }
        
        std::cerr << "  高电平切线:\n";
        for(const auto& vline : plan.vertical_lines) {
          if(vline.is_high) {
            std::cerr << "    竖线: x=" << vline.coord/2 << ", y∈[" 
                      << vline.span_min2/2 << "," << vline.span_max2/2 << "]\n";
          }
        }
        for(const auto& hline : plan.horizontal_lines) {
          if(hline.is_high) {
            std::cerr << "    横线: y=" << hline.coord/2 << ", x∈[" 
                      << hline.span_min2/2 << "," << hline.span_max2/2 << "]\n";
          }
        }
        
        debug_count++;
        break;  // 只看这个AA的第一个孤立片段
      }
    }
  } else {
    std::cerr << "\n✓ 所有孤立片段都是合理的（只有低电平切线或无切线）\n";
  }
  
  // 统计高电平/低电平PO的bbox分布
  int high_po_count = 0, low_po_count = 0;
  int high_po_minx = INT32_MAX, high_po_miny = INT32_MAX, high_po_maxx = INT32_MIN, high_po_maxy = INT32_MIN;
  int low_po_minx = INT32_MAX, low_po_miny = INT32_MAX, low_po_maxx = INT32_MIN, low_po_maxy = INT32_MIN;
  
  if(gate_ctx.poly_lid != UINT32_MAX) {
    const auto& polys = L->layers[gate_ctx.poly_lid].polys;
    for(u32 pid = 0; pid < polys.size(); pid++) {
      const auto& poly = polys[pid];
      bool is_high = (pid < gate_ctx.poly_high.size() && gate_ctx.poly_high[pid]);
      
      if(is_high) {
        high_po_count++;
        high_po_minx = std::min(high_po_minx, poly.bb.minx);
        high_po_miny = std::min(high_po_miny, poly.bb.miny);
        high_po_maxx = std::max(high_po_maxx, poly.bb.maxx);
        high_po_maxy = std::max(high_po_maxy, poly.bb.maxy);
      } else {
        low_po_count++;
        low_po_minx = std::min(low_po_minx, poly.bb.minx);
        low_po_miny = std::min(low_po_miny, poly.bb.miny);
        low_po_maxx = std::max(low_po_maxx, poly.bb.maxx);
        low_po_maxy = std::max(low_po_maxy, poly.bb.maxy);
      }
    }
  }
  
  std::cerr << "\n=== 高/低电平PO空间分布 ===\n";
  std::cerr << "高电平PO: " << high_po_count << " 个\n";
  if(high_po_count > 0) {
    std::cerr << "  范围: [" << high_po_minx << "," << high_po_miny << "] 到 [" 
              << high_po_maxx << "," << high_po_maxy << "]\n";
    std::cerr << "  宽度=" << (high_po_maxx - high_po_minx) << ", 高度=" << (high_po_maxy - high_po_miny) << "\n";
  }
  std::cerr << "低电平PO: " << low_po_count << " 个\n";
  if(low_po_count > 0) {
    std::cerr << "  范围: [" << low_po_minx << "," << low_po_miny << "] 到 [" 
              << low_po_maxx << "," << low_po_maxy << "]\n";
    std::cerr << "  宽度=" << (low_po_maxx - low_po_minx) << ", 高度=" << (low_po_maxy - low_po_miny) << "\n";
  }
  
  // 统计Phase 2访问的AA的bbox分布
  int visited_aa_count = 0;
  int visited_aa_minx = INT32_MAX, visited_aa_miny = INT32_MAX, visited_aa_maxx = INT32_MIN, visited_aa_maxy = INT32_MIN;
  
  if(gate_ctx.aa_lid != UINT32_MAX) {
    const auto& aa_polys = L->layers[gate_ctx.aa_lid].polys;
    for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
      if(slices.ready && aa_pid < aa_polys.size()) {
        visited_aa_count++;
        const auto& aa = aa_polys[aa_pid];
        visited_aa_minx = std::min(visited_aa_minx, aa.bb.minx);
        visited_aa_miny = std::min(visited_aa_miny, aa.bb.miny);
        visited_aa_maxx = std::max(visited_aa_maxx, aa.bb.maxx);
        visited_aa_maxy = std::max(visited_aa_maxy, aa.bb.maxy);
      }
    }
  }
  
  std::cerr << "\n=== Phase 2访问的AA空间分布 ===\n";
  std::cerr << "访问的AA: " << visited_aa_count << " 个（总AA数: " << L->layers[gate_ctx.aa_lid].polys.size() << "）\n";
  if(visited_aa_count > 0) {
    std::cerr << "  范围: [" << visited_aa_minx << "," << visited_aa_miny << "] 到 [" 
              << visited_aa_maxx << "," << visited_aa_maxy << "]\n";
    std::cerr << "  宽度=" << (visited_aa_maxx - visited_aa_minx) << ", 高度=" << (visited_aa_maxy - visited_aa_miny) << "\n";
  }
  
  // 打印s1和s2的位置（从规则中读取）
  std::cerr << "\n=== 起点位置 ===\n";
  std::cerr << "s1: layer_id=" << lid1 << ", pid=" << pid1;
  if(lid1 < L->layers.size() && pid1 < L->layers[lid1].polys.size()) {
    const auto& s1_bb = L->layers[lid1].polys[pid1].bb;
    std::cerr << ", bbox中心=(" << (s1_bb.minx + s1_bb.maxx)/2 << "," << (s1_bb.miny + s1_bb.maxy)/2 << ")";
  }
  std::cerr << "\n";
  std::cerr << "s2: layer_id=" << lid2 << ", pid=" << pid2;
  if(lid2 < L->layers.size() && pid2 < L->layers[lid2].polys.size()) {
    const auto& s2_bb = L->layers[lid2].polys[pid2].bb;
    std::cerr << ", bbox中心=(" << (s2_bb.minx + s2_bb.maxx)/2 << "," << (s2_bb.miny + s2_bb.maxy)/2 << ")";
  }
  std::cerr << "\n";
  
  return out;
}

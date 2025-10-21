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
  {
    std::vector<u32> cand; cand.reserve(128);
    grids[lid].query(U.bb, cand);
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    for(u32 vpid : cand){
      if(vpid == U.pid) continue;
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
          // Phase 1：实验-回退改进1，当作普通多边形处理
          // 允许通过AA扩张到其他层
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
  return out;
}

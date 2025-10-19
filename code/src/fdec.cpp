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
  for(const auto& chain : R->vias){
    for(size_t i=1;i<chain.size();++i){
      auto it1=n2id.find(chain[i-1]); auto it2=n2id.find(chain[i]);
      if(it1!=n2id.end() && it2!=n2id.end()) add_adj(it1->second, it2->second);
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
                                   std::unordered_set<u64>& visited) const{
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
      
      // 特殊处理：如果nlid是AA层，使用入口片逻辑
      if(nlid == gate_ctx.aa_lid && gate_ctx.aa_lid != UINT32_MAX) {
        // 先切割AA（如果还没切）
        const_cast<GateCtx&>(gate_ctx).execute_cut(*L, npid);
        // 只入队入口片，并在AA内沿高切线泛洪
        gate_ctx.enqueue_entry_slices_from(gid, npid, q, visited, *L);
        continue; // 不再把整块AA入队
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
    expand_frontier_no_gate(u, q, visited);
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

  // 阶段1：s1 BFS（访问到 poly_layer 即置 high）
  {
    std::unordered_set<u64> visited1; visited1.reserve(1<<20);
    std::queue<u64> q1;
    u64 g1 = make_gid(lid1,pid1);
    visited1.insert(g1); q1.push(g1);
    
    int poly_high_count = 0;

    while(!q1.empty()){
      u64 u=q1.front(); q1.pop();
      u32 lid=gid_lid(u), pid=gid_pid(u);
      if(lid==gate_ctx.poly_lid) {
        gate_ctx.poly_high[pid]=1;
        poly_high_count++;
      }
      // 同"无 gate"扩张：同层 + via 相邻层
      expand_frontier_no_gate(u, q1, visited1);
    }
    
    std::cerr << "Phase 1: visited " << visited1.size() << " polygons, marked " << poly_high_count << " poly as high\n";
  }

  // 阶段2：s2 BFS（命中 AA 时执行多次切割，并通过导通边扩展）
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

  u64 g2 = make_gid(lid2,pid2);
  visited2.insert(g2); q2.push(g2);

  while(!q2.empty()){
    u64 u=q2.front(); q2.pop();
    
    // 处理切割后的片段
    if(gate_ctx.is_slice_gid(u)){
      out.push_back(u);
      total_pieces++;
      
      // 常规扩张（到其他层，AA的内部扩张已在enqueue_entry_slices_from中完成）
      expand_frontier_no_gate(u, q2, visited2);
      continue;
    }

    // 统计切割数量
    u32 lid = gid_lid(u), pid = gid_pid(u);
    if(lid == gate_ctx.aa_lid) {
      auto it_slices = gate_ctx.aa_slices.find(pid);
      if(it_slices != gate_ctx.aa_slices.end() && it_slices->second.ready) {
        if(it_slices->second.pieces.size() > 1) {
          aa_cut_count++;
        }
      }
    }

    // 其他层的正常处理
    out.push_back(u);
    expand_frontier_no_gate(u, q2, visited2);
  }
  
  // 统计有多少片段没被访问
  int total_generated_pieces = 0;
  for(const auto& [aa_pid, slices] : gate_ctx.aa_slices) {
    if(slices.ready) {
      total_generated_pieces += slices.pieces.size();
    }
  }
  
  std::cerr << "Phase 2: cut " << aa_cut_count << " AA polygons\n";
  std::cerr << "  Total generated pieces: " << total_generated_pieces << "\n";
  std::cerr << "  Total visited pieces: " << total_pieces << "\n";
  std::cerr << "  Unvisited pieces: " << (total_generated_pieces - total_pieces) << "\n";
  return out;
}

#pragma once
#include "types.hpp"
#include "grid.hpp"
#include "gate.hpp"
#include <queue>
#include <unordered_set>

struct FDEC {
  const Layout* L=nullptr;
  const Rule*   R=nullptr;
  Config        cfg;

  // one grid per layer
  std::vector<Grid> grids;

  // via adjacency: lid -> list of adjacent lids
  std::vector<std::vector<u32>> adj_layers;

  mutable GateCtx gate_ctx; // 缓存 Gate 阶段运行时数据（Lazy-Cut 切片等）

  // visited：先用哈希集合替代位图，便于骨架可跑（后续可替换为位图/roaring）
  // 为简洁放在实现文件里使用 unordered_set<u64>

  void build(const Layout& layout, const Rule& rule, const Config& cfg);

  // find poly pid that contains seed point on specified layer (boundary counts as inside)
  u32 locate_seed_pid(u32 lid, const Pt& p) const;

  // First/Second task (no gate): multi-source BFS
  std::vector<u64> trace_no_gate(const std::vector<std::pair<u32,u32>>& seed_lid_pid) const;

  // Third task (with gate, two-phase: s1 -> poly_high; s2 -> BFS+lazy-cut)
  std::vector<u64> trace_with_gate() const;

  const GateCtx* gate_result() const { return gate_ctx.poly_lid==UINT32_MAX ? nullptr : &gate_ctx; }

private:
  void expand_frontier_no_gate(u64 gid, std::queue<u64>& q,
                               std::unordered_set<u64>& visited) const;
  void parse_via_to_adjacency();
};

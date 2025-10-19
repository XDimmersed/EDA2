#pragma once
#include "types.hpp"
#include "grid.hpp"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <queue>

struct GateCtx {
  u32 poly_lid=UINT32_MAX, aa_lid=UINT32_MAX;
  std::vector<uint8_t> poly_high;                 // poly pid -> 0/1
  std::vector<std::vector<u32>> poly_to_aa_cands; // coarse candidates per poly
  std::vector<std::vector<u32>> aa_to_poly_cands; // coarse candidates per AA
  Grid poly_grid, aa_grid;                        // grids for coarse query reuse
  
  // 切割线定义
  struct CutLine {
    i32 coord;      // x坐标（竖切）或y坐标（横切），使用2×坐标
    bool is_high;   // 高电平（导通）还是低电平（断开）
    i32 span_min2;  // 作用区间（另一维度）的下界，使用2×坐标
    i32 span_max2;  // 作用区间（另一维度）的上界，使用2×坐标
  };
  
  // AA切割计划
  struct AAPlan {
    std::vector<CutLine> vertical_lines;   // 竖切线（按x排序）
    std::vector<CutLine> horizontal_lines; // 横切线（按y排序）
    bool computed = false;
  };
  
  // AA切割后的多个片段
  struct AASlices {
    std::vector<Poly> pieces;           // 所有片段
    std::vector<std::vector<u32>> adjacency; // 导通边：pieces[i]的邻居索引列表
    bool ready = false;
  };
  
  std::unordered_map<u32, AAPlan> aa_plans;     // aa pid -> 切割计划
  std::unordered_map<u32, AASlices> aa_slices;  // aa pid -> 切割后的片段

  void init_empty(){
    poly_lid=aa_lid=UINT32_MAX;
    poly_high.clear();
    poly_to_aa_cands.clear();
    aa_to_poly_cands.clear();
    aa_plans.clear();
    aa_slices.clear();
    poly_grid = Grid();
    aa_grid = Grid();
  }

  // Build coarse candidates (bbox + grid probe)
  void build_candidates(const Layout& L, int cell_hint);
  
  // 计算AA的切割计划（收集所有贯穿的Poly）
  void compute_aa_plan(const Layout& L, u32 aa_pid);
  
  // 执行切割并建立导通边
  void execute_cut(const Layout& L, u32 aa_pid);
  
  // 从U扩展到AA时，只入队入口片，并在AA内沿高切线泛洪
  void enqueue_entry_slices_from(u64 u_gid, u32 aa_pid, 
                                  std::queue<u64>& q, 
                                  std::unordered_set<u64>& vis,
                                  const Layout& L) const;
  
  // 获取片段的gid（使用注册表避免8位限制）
  struct SliceKey {
    u32 aa_pid;
    u32 piece_idx;
    bool operator==(const SliceKey& o) const { 
      return aa_pid == o.aa_pid && piece_idx == o.piece_idx; 
    }
  };
  
  struct SliceKeyHash {
    size_t operator()(const SliceKey& k) const {
      return ((size_t)k.aa_pid << 32) | k.piece_idx;
    }
  };
  
  mutable std::unordered_map<SliceKey, u64, SliceKeyHash> slice_key_to_gid;
  mutable std::unordered_map<u64, SliceKey> gid_to_slice_key;
  mutable u64 next_slice_gid = (1ULL << 62); // 使用高位避免与普通gid冲突
  
  u64 make_slice_gid(u32 aa_pid, u32 piece_idx) const {
    SliceKey key{aa_pid, piece_idx};
    auto it = slice_key_to_gid.find(key);
    if(it != slice_key_to_gid.end()) return it->second;
    
    u64 gid = next_slice_gid++;
    slice_key_to_gid[key] = gid;
    gid_to_slice_key[gid] = key;
    return gid;
  }
  
  SliceKey gid_to_key(u64 gid) const {
    auto it = gid_to_slice_key.find(gid);
    if(it != gid_to_slice_key.end()) return it->second;
    return {UINT32_MAX, UINT32_MAX};
  }
  
  bool is_slice_gid(u64 gid) const { 
    return gid >= (1ULL << 62); // 切片gid在高位区域
  }
  
  const Poly* slice_from_gid(u64 gid) const {
    auto key = gid_to_key(gid);
    if(key.aa_pid == UINT32_MAX) return nullptr;
    auto it = aa_slices.find(key.aa_pid);
    if(it == aa_slices.end() || !it->second.ready) return nullptr;
    if(key.piece_idx >= it->second.pieces.size()) return nullptr;
    return &it->second.pieces[key.piece_idx];
  }
};

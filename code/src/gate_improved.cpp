#include "gate.hpp"
#include "geom.hpp"

#include <algorithm>
#include <utility>
#include <iostream>
#include <map>

namespace {

BBox bbox_intersection(const BBox& a, const BBox& b){
  BBox r;
  r.minx = std::max(a.minx, b.minx);
  r.maxx = std::min(a.maxx, b.maxx);
  r.miny = std::max(a.miny, b.miny);
  r.maxy = std::min(a.maxy, b.maxy);
  if(r.minx > r.maxx || r.miny > r.maxy){
    r.minx=r.maxx=r.miny=r.maxy=0;
  }
  return r;
}

// 改进的贯穿判据：检查poly是否"穿过"AA
bool is_full_penetration(const Poly& poly, const Poly& aa) {
  BBox intersection = bbox_intersection(poly.bb, aa.bb);
  
  // 如果交集为空，肯定不贯穿
  if(intersection.minx >= intersection.maxx || intersection.miny >= intersection.maxy) {
    return false;
  }
  
  // 计算交集的面积和AA的面积
  long long intersection_area = (long long)(intersection.maxx - intersection.minx) * 
                               (long long)(intersection.maxy - intersection.miny);
  long long aa_area = (long long)(aa.bb.maxx - aa.bb.minx) * 
                     (long long)(aa.bb.maxy - aa.bb.miny);
  
  // 如果交集面积占AA面积的很大比例（比如>80%），认为完全贯穿
  if(aa_area > 0 && intersection_area * 100 > aa_area * 80) {
    return true;
  }
  
  // 检查是否真正贯穿
  // 横向贯穿：交集必须接触AA的上下边界
  bool horizontal_penetration = (intersection.miny <= aa.bb.miny && 
                                intersection.maxy >= aa.bb.maxy);
  
  // 纵向贯穿：交集必须接触AA的左右边界
  bool vertical_penetration = (intersection.minx <= aa.bb.minx && 
                              intersection.maxx >= aa.bb.maxx);
  
  // 放宽条件：如果交集在某个方向上几乎覆盖AA，也认为贯穿
  if(horizontal_penetration) {
    // 检查poly是否在AA的左右两侧都有延伸（允许小误差）
    bool extends_left = (poly.bb.minx < aa.bb.minx + 100);  // 允许100单位的误差
    bool extends_right = (poly.bb.maxx > aa.bb.maxx - 100);
    if(extends_left && extends_right) return true;
  }
  
  if(vertical_penetration) {
    // 检查poly是否在AA的上下两侧都有延伸（允许小误差）
    bool extends_bottom = (poly.bb.miny < aa.bb.miny + 100);  // 允许100单位的误差
    bool extends_top = (poly.bb.maxy > aa.bb.maxy - 100);
    if(extends_bottom && extends_top) return true;
  }
  
  return horizontal_penetration || vertical_penetration;
}

// 改进的Sutherland-Hodgman裁剪算法，支持半格切线
struct ClipResult {
  std::vector<Pt> vertices;
  bool is_valid;
  
  ClipResult() : is_valid(false) {}
  ClipResult(std::vector<Pt> v) : vertices(std::move(v)), is_valid(vertices.size() >= 3) {}
};

// 完整的Sutherland-Hodgman裁剪
ClipResult sutherland_hodgman_clip(const std::vector<Pt>& subject, 
                                   bool is_horizontal, 
                                   i32 coord2, 
                                   bool keep_lower) {
  if(subject.size() < 3) return ClipResult();
  
  std::vector<Pt> result = subject;
  
  // 定义裁剪平面
  auto inside = [&](const Pt& p) -> bool {
    if(is_horizontal) {
      return keep_lower ? (p.y * 2 < coord2) : (p.y * 2 >= coord2);
    } else {
      return keep_lower ? (p.x * 2 < coord2) : (p.x * 2 >= coord2);
    }
  };
  
  auto intersect = [&](const Pt& p1, const Pt& p2) -> Pt {
    if(is_horizontal) {
      // 横切：计算与y=coord2/2的交点
      if(p1.y == p2.y) return p1; // 水平边
      i32 y_cut = coord2 / 2;
      if(coord2 % 2 == 0) {
        // 整数切线
        long long x = p1.x + (long long)(y_cut - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
        return {(i32)x, y_cut};
      } else {
        // 半格切线，返回边的中点
        return {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2};
      }
    } else {
      // 竖切：计算与x=coord2/2的交点
      if(p1.x == p2.x) return p1; // 垂直边
      i32 x_cut = coord2 / 2;
      if(coord2 % 2 == 0) {
        // 整数切线
        long long y = p1.y + (long long)(x_cut - p1.x) * (p2.y - p1.y) / (p2.x - p1.x);
        return {x_cut, (i32)y};
      } else {
        // 半格切线，返回边的中点
        return {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2};
      }
    }
  };
  
  std::vector<Pt> input = result;
  result.clear();
  
  if(input.empty()) return ClipResult();
  
  Pt s = input.back();
  for(const Pt& e : input) {
    if(inside(e)) {
      if(!inside(s)) {
        result.push_back(intersect(s, e));
      }
      result.push_back(e);
    } else if(inside(s)) {
      result.push_back(intersect(s, e));
    }
    s = e;
  }
  
  // 去重相邻的重复点
  if(!result.empty()) {
    std::vector<Pt> deduped;
    deduped.push_back(result[0]);
    for(size_t i = 1; i < result.size(); ++i) {
      if(result[i].x != result[i-1].x || result[i].y != result[i-1].y) {
        deduped.push_back(result[i]);
      }
    }
    result = std::move(deduped);
  }
  
  return ClipResult(result);
}

// 规范化多边形：确保CCW，重建H/V边数组，更新bbox
void normalize_polygon(Poly& poly) {
  if(poly.v.size() < 3) return;
  
  // 确保逆时针
  ensure_ccw(poly.v);
  
  // 重建H/V边数组
  poly.H.clear();
  poly.V.clear();
  
  for(size_t i = 0; i < poly.v.size(); ++i) {
    const Pt& p1 = poly.v[i];
    const Pt& p2 = poly.v[(i + 1) % poly.v.size()];
    
    if(p1.y == p2.y) {
      // 水平边
      poly.H.push_back({p1, p2});
    } else if(p1.x == p2.x) {
      // 垂直边
      poly.V.push_back({p1, p2});
    }
  }
  
  // 更新bbox
  poly.bb = bbox_of(poly.v);
}

// 切割结果结构，包含左右/上下片段和导通信息
struct CutResult {
  std::vector<Poly> pieces;
  std::vector<std::vector<u32>> adjacency; // 导通边
  bool is_high; // 是否高电平导通
  
  CutResult() : is_high(false) {}
};

// 改进的切割函数，在切割时就建立导通关系
CutResult cut_polygon_with_adjacency(const Poly& original, 
                                     bool is_horizontal, 
                                     i32 coord2, 
                                     bool is_high) {
  CutResult result;
  result.is_high = is_high;
  
  if(original.v.size() < 3) return result;
  
  // 使用改进的Sutherland-Hodgman裁剪
  auto left_result = sutherland_hodgman_clip(original.v, is_horizontal, coord2, true);
  auto right_result = sutherland_hodgman_clip(original.v, is_horizontal, coord2, false);
  
  if(left_result.is_valid) {
    Poly left_piece;
    left_piece.v = std::move(left_result.vertices);
    left_piece.lid = original.lid;
    normalize_polygon(left_piece);
    result.pieces.push_back(std::move(left_piece));
  }
  
  if(right_result.is_valid) {
    Poly right_piece;
    right_piece.v = std::move(right_result.vertices);
    right_piece.lid = original.lid;
    normalize_polygon(right_piece);
    result.pieces.push_back(std::move(right_piece));
  }
  
  // 如果高电平且有两个片段，建立导通边
  if(is_high && result.pieces.size() == 2) {
    result.adjacency.resize(2);
    result.adjacency[0].push_back(1);
    result.adjacency[1].push_back(0);
  }
  
  return result;
}

} // namespace

void GateCtx::build_candidates(const Layout& L, int cell_hint){
  if(poly_lid >= L.layers.size() || aa_lid >= L.layers.size()) return;
  
  const auto& poly_layer = L.layers[poly_lid];
  const auto& aa_layer = L.layers[aa_lid];
  
  poly_high.assign(poly_layer.polys.size(), 0);
  poly_to_aa_cands.assign(poly_layer.polys.size(), {});
  aa_to_poly_cands.assign(aa_layer.polys.size(), {});
  
  // 构建网格
  poly_grid = Grid::build(poly_layer.polys, cell_hint);
  aa_grid = Grid::build(aa_layer.polys, cell_hint);
  
  // 建立候选关系
  for(u32 poly_pid = 0; poly_pid < poly_layer.polys.size(); ++poly_pid) {
    const auto& poly = poly_layer.polys[poly_pid];
    std::vector<u32> candidates;
    aa_grid.query(poly.bb, candidates);
    
    for(u32 aa_pid : candidates) {
      const auto& aa = aa_layer.polys[aa_pid];
      if(bbox_overlap(poly.bb, aa.bb)) {
        poly_to_aa_cands[poly_pid].push_back(aa_pid);
        aa_to_poly_cands[aa_pid].push_back(poly_pid);
      }
    }
  }
}

void GateCtx::compute_aa_plan(const Layout& L, u32 aa_pid) {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  
  auto& plan = aa_plans[aa_pid];
  if(plan.computed) return;
  
  const auto& aa = L.layers[aa_lid].polys[aa_pid];
  plan.vertical_lines.clear();
  plan.horizontal_lines.clear();
  
  // 收集所有高电平Poly的贯穿线
  for(u32 poly_pid : aa_to_poly_cands[aa_pid]) {
    if(poly_pid >= poly_high.size() || !poly_high[poly_pid]) continue;
    
    const auto& poly = L.layers[poly_lid].polys[poly_pid];
    if(!is_full_penetration(poly, aa)) continue;
    
    BBox intersection = bbox_intersection(poly.bb, aa.bb);
    
    // 检查横向贯穿
    if(intersection.miny <= aa.bb.miny && intersection.maxy >= aa.bb.maxy) {
      // 使用中线作为切割线
      i32 y_cut = (intersection.miny + intersection.maxy) / 2;
      plan.horizontal_lines.push_back({y_cut * 2, true}); // 2×坐标
    }
    
    // 检查纵向贯穿
    if(intersection.minx <= aa.bb.minx && intersection.maxx >= aa.bb.maxx) {
      // 使用中线作为切割线
      i32 x_cut = (intersection.minx + intersection.maxx) / 2;
      plan.vertical_lines.push_back({x_cut * 2, true}); // 2×坐标
    }
  }
  
  // 排序去重
  std::sort(plan.vertical_lines.begin(), plan.vertical_lines.end(), 
            [](const CutLine& a, const CutLine& b) { return a.coord < b.coord; });
  std::sort(plan.horizontal_lines.begin(), plan.horizontal_lines.end(), 
            [](const CutLine& a, const CutLine& b) { return a.coord < b.coord; });
  
  plan.computed = true;
}

void GateCtx::execute_cut(const Layout& L, u32 aa_pid) {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  
  auto& slices = aa_slices[aa_pid];
  if(slices.ready) return;
  
  compute_aa_plan(L, aa_pid);
  const auto& plan = aa_plans[aa_pid];
  
  if(plan.vertical_lines.empty() && plan.horizontal_lines.empty()) {
    // 没有切割线，保持原样
    slices.pieces.push_back(L.layers[aa_lid].polys[aa_pid]);
    slices.adjacency.resize(1);
    slices.ready = true;
    return;
  }
  
  // 开始切割
  std::vector<Poly> current_pieces = {L.layers[aa_lid].polys[aa_pid]};
  std::vector<std::vector<u32>> current_adjacency = {{}};
  
  // 先竖切
  for(const auto& line : plan.vertical_lines) {
    std::vector<Poly> new_pieces;
    std::vector<std::vector<u32>> new_adjacency;
    
    for(size_t i = 0; i < current_pieces.size(); ++i) {
      const auto& piece = current_pieces[i];
      
      if(piece.bb.minx * 2 < line.coord && line.coord < piece.bb.maxx * 2) {
        // 需要切割
        auto cut_result = cut_polygon_with_adjacency(piece, false, line.coord, line.is_high);
        
        if(cut_result.pieces.size() == 2) {
          // 成功切割成两片
          size_t left_idx = new_pieces.size();
          size_t right_idx = new_pieces.size() + 1;
          
          new_pieces.push_back(std::move(cut_result.pieces[0]));
          new_pieces.push_back(std::move(cut_result.pieces[1]));
          
          // 建立导通边
          new_adjacency.resize(new_pieces.size());
          if(cut_result.is_high) {
            new_adjacency[left_idx].push_back(right_idx);
            new_adjacency[right_idx].push_back(left_idx);
          }
          
          // 继承原片段的邻接关系
          for(u32 neighbor : current_adjacency[i]) {
            new_adjacency[left_idx].push_back(neighbor);
            new_adjacency[neighbor].push_back(left_idx);
            new_adjacency[right_idx].push_back(neighbor);
            new_adjacency[neighbor].push_back(right_idx);
          }
        } else {
          // 切割失败，保持原样
          new_pieces.push_back(piece);
          new_adjacency.resize(new_pieces.size());
          new_adjacency[new_pieces.size()-1] = current_adjacency[i];
        }
      } else {
        // 不需要切割
        new_pieces.push_back(piece);
        new_adjacency.resize(new_pieces.size());
        new_adjacency[new_pieces.size()-1] = current_adjacency[i];
      }
    }
    
    current_pieces = std::move(new_pieces);
    current_adjacency = std::move(new_adjacency);
  }
  
  // 再横切
  for(const auto& line : plan.horizontal_lines) {
    std::vector<Poly> new_pieces;
    std::vector<std::vector<u32>> new_adjacency;
    
    for(size_t i = 0; i < current_pieces.size(); ++i) {
      const auto& piece = current_pieces[i];
      
      if(piece.bb.miny * 2 < line.coord && line.coord < piece.bb.maxy * 2) {
        // 需要切割
        auto cut_result = cut_polygon_with_adjacency(piece, true, line.coord, line.is_high);
        
        if(cut_result.pieces.size() == 2) {
          // 成功切割成两片
          size_t lower_idx = new_pieces.size();
          size_t upper_idx = new_pieces.size() + 1;
          
          new_pieces.push_back(std::move(cut_result.pieces[0]));
          new_pieces.push_back(std::move(cut_result.pieces[1]));
          
          // 建立导通边
          new_adjacency.resize(new_pieces.size());
          if(cut_result.is_high) {
            new_adjacency[lower_idx].push_back(upper_idx);
            new_adjacency[upper_idx].push_back(lower_idx);
          }
          
          // 继承原片段的邻接关系
          for(u32 neighbor : current_adjacency[i]) {
            new_adjacency[lower_idx].push_back(neighbor);
            new_adjacency[neighbor].push_back(lower_idx);
            new_adjacency[upper_idx].push_back(neighbor);
            new_adjacency[neighbor].push_back(upper_idx);
          }
        } else {
          // 切割失败，保持原样
          new_pieces.push_back(piece);
          new_adjacency.resize(new_pieces.size());
          new_adjacency[new_pieces.size()-1] = current_adjacency[i];
        }
      } else {
        // 不需要切割
        new_pieces.push_back(piece);
        new_adjacency.resize(new_pieces.size());
        new_adjacency[new_pieces.size()-1] = current_adjacency[i];
      }
    }
    
    current_pieces = std::move(new_pieces);
    current_adjacency = std::move(new_adjacency);
  }
  
  // 保存结果
  slices.pieces = std::move(current_pieces);
  slices.adjacency = std::move(current_adjacency);
  slices.ready = true;
}

void GateCtx::enqueue_entry_slices_from(u64 u_gid, u32 aa_pid, 
                                        std::queue<u64>& q, 
                                        std::unordered_set<u64>& vis,
                                        const Layout& L) const {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  
  auto it = aa_slices.find(aa_pid);
  if(it == aa_slices.end() || !it->second.ready) return;
  
  const auto& slices = it->second;
  const auto& u_poly = L.layers[gid_lid(u_gid)].polys[gid_pid(u_gid)];
  
  // 找到与u_poly相交的片段
  for(size_t i = 0; i < slices.pieces.size(); ++i) {
    const auto& slice = slices.pieces[i];
    if(bbox_overlap(u_poly.bb, slice.bb) && poly_intersect_manhattan(u_poly, slice)) {
      u64 slice_gid = make_slice_gid(aa_pid, i);
      if(vis.find(slice_gid) == vis.end()) {
        vis.insert(slice_gid);
        q.push(slice_gid);
        
        // 沿导通边扩展
        for(u32 neighbor_idx : slices.adjacency[i]) {
          u64 neighbor_gid = make_slice_gid(aa_pid, neighbor_idx);
          if(vis.find(neighbor_gid) == vis.end()) {
            vis.insert(neighbor_gid);
            q.push(neighbor_gid);
          }
        }
      }
    }
  }
}
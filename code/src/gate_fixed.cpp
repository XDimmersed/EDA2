#include "gate.hpp"
#include "geom.hpp"
#include "parser.hpp"
#include <algorithm>
#include <map>
#include <queue>
#include <unordered_set>

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

// 真正的贯穿判据：检查Poly是否真正在AA内形成通道
bool is_full_penetration(const Poly& poly, const Poly& aa) {
  // 首先检查bbox交集
  BBox intersection = bbox_intersection(poly.bb, aa.bb);
  if(intersection.minx >= intersection.maxx || intersection.miny >= intersection.maxy) {
    return false;
  }
  
  // 计算交集面积
  long long intersection_area = (long long)(intersection.maxx - intersection.minx) * 
                               (long long)(intersection.maxy - intersection.miny);
  long long aa_area = (long long)(aa.bb.maxx - aa.bb.minx) * 
                     (long long)(aa.bb.maxy - aa.bb.miny);
  
  // 如果交集面积占AA面积的很大比例（>85%），认为完全贯穿
  if(aa_area > 0 && intersection_area * 100 > aa_area * 85) {
    return true;
  }
  
  // 真正的贯穿检查：Poly必须真正跨越AA的边界
  // 横向贯穿：Poly必须与AA的上下边界都有交点，且在AA内部形成连通区域
  bool horizontal_penetration = false;
  if(poly.bb.miny <= aa.bb.miny && poly.bb.maxy >= aa.bb.maxy) {
    // Poly的y范围包含AA的y范围，检查是否在x方向上有延伸
    if(poly.bb.minx < aa.bb.minx && poly.bb.maxx > aa.bb.maxx) {
      horizontal_penetration = true;
    }
  }
  
  // 纵向贯穿：Poly必须与AA的左右边界都有交点，且在AA内部形成连通区域
  bool vertical_penetration = false;
  if(poly.bb.minx <= aa.bb.minx && poly.bb.maxx >= aa.bb.maxx) {
    // Poly的x范围包含AA的x范围，检查是否在y方向上有延伸
    if(poly.bb.miny < aa.bb.miny && poly.bb.maxy > aa.bb.maxy) {
      vertical_penetration = true;
    }
  }
  
  return horizontal_penetration || vertical_penetration;
}

// 完整的Sutherland-Hodgman裁剪算法，正确处理半格坐标
struct ClipResult {
  std::vector<Pt> vertices;
  bool success;
  
  ClipResult() : success(false) {}
  ClipResult(const std::vector<Pt>& v) : vertices(v), success(v.size() >= 3) {}
};

ClipResult sutherland_hodgman_clip(const std::vector<Pt>& subject,
                                   bool is_horizontal,
                                   i32 coord2,
                                   bool keep_lower) {
  if(subject.size() < 3) return ClipResult();
  
  // 定义裁剪平面（在2×坐标空间）
  auto inside = [&](const Pt& p) -> bool {
    if(is_horizontal) {
      return keep_lower ? (p.y * 2 < coord2) : (p.y * 2 >= coord2);
    } else {
      return keep_lower ? (p.x * 2 < coord2) : (p.x * 2 >= coord2);
    }
  };
  
  // 计算精确的交点
  auto intersect = [&](const Pt& p1, const Pt& p2) -> Pt {
    if(is_horizontal) {
      // 横切：计算与y=coord2/2的交点
      if(p1.y == p2.y) return p1; // 水平边
      
      // 在2×坐标空间计算交点
      i32 y1_2 = p1.y * 2, y2_2 = p2.y * 2;
      if(y1_2 == y2_2) return p1; // 水平边
      
      // 计算交点：x = x1 + (coord2 - y1_2) * (x2 - x1) / (y2_2 - y1_2)
      long long x = p1.x + (long long)(coord2 - y1_2) * (p2.x - p1.x) / (y2_2 - y1_2);
      i32 y = coord2 / 2; // 精确的y坐标
      
      return {(i32)x, y};
    } else {
      // 竖切：计算与x=coord2/2的交点
      if(p1.x == p2.x) return p1; // 垂直边
      
      // 在2×坐标空间计算交点
      i32 x1_2 = p1.x * 2, x2_2 = p2.x * 2;
      if(x1_2 == x2_2) return p1; // 垂直边
      
      // 计算交点：y = y1 + (coord2 - x1_2) * (y2 - y1) / (x2_2 - x1_2)
      long long y = p1.y + (long long)(coord2 - x1_2) * (p2.y - p1.y) / (x2_2 - x1_2);
      i32 x = coord2 / 2; // 精确的x坐标
      
      return {x, (i32)y};
    }
  };
  
  std::vector<Pt> result;
  
  if(subject.empty()) return ClipResult();
  
  Pt s = subject.back();
  for(const Pt& e : subject) {
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
  
  // 强制闭环：确保首尾点相同
  if(!result.empty() && (result[0].x != result.back().x || result[0].y != result.back().y)) {
    result.push_back(result[0]);
  }
  
  return ClipResult(result);
}

// 规范化多边形：确保CCW，重建H/V边数组，更新bbox
void normalize_polygon(Poly& poly) {
  if(poly.v.size() < 3) return;
  
  // 确保逆时针
  ensure_ccw(poly.v);
  
  // 重建H和V边数组
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

// 切割多边形并建立邻接关系
struct CutResult {
  std::vector<Poly> pieces;
  bool is_high;
  
  CutResult() : is_high(false) {}
  CutResult(const std::vector<Poly>& p, bool high) : pieces(p), is_high(high) {}
};

CutResult cut_polygon_with_adjacency(const Poly& original,
                                     bool is_horizontal,
                                     i32 coord2,
                                     bool is_high) {
  // 使用完整的Sutherland-Hodgman算法进行切割
  auto clip_result = sutherland_hodgman_clip(original.v, is_horizontal, coord2, true);
  if(!clip_result.success) {
    return CutResult();
  }
  
  // 创建左/下片段
  Poly left_piece;
  left_piece.v = clip_result.vertices;
  left_piece.lid = original.lid;
  left_piece.pid = original.pid;
  left_piece.gid = original.gid;
  normalize_polygon(left_piece);
  
  // 创建右/上片段
  auto right_clip = sutherland_hodgman_clip(original.v, is_horizontal, coord2, false);
  if(!right_clip.success) {
    return CutResult();
  }
  
  Poly right_piece;
  right_piece.v = right_clip.vertices;
  right_piece.lid = original.lid;
  right_piece.pid = original.pid;
  right_piece.gid = original.gid;
  normalize_polygon(right_piece);
  
  return CutResult({left_piece, right_piece}, is_high);
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
  
  const Poly& aa = L.layers[aa_lid].polys[aa_pid];
  
  // 使用map来合并相同坐标的切割线，高优先级覆盖低优先级
  std::map<i32, bool> vertical_cuts, horizontal_cuts;
  
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
      horizontal_cuts[y_cut * 2] = true; // 高电平
    }
    
    // 检查纵向贯穿
    if(intersection.minx <= aa.bb.minx && intersection.maxx >= aa.bb.maxx) {
      // 使用中线作为切割线
      i32 x_cut = (intersection.minx + intersection.maxx) / 2;
      vertical_cuts[x_cut * 2] = true; // 高电平
    }
  }
  
  // 转换为排序的向量
  for(const auto& [coord, is_high] : vertical_cuts) {
    plan.vertical_lines.push_back({coord, is_high});
  }
  for(const auto& [coord, is_high] : horizontal_cuts) {
    plan.horizontal_lines.push_back({coord, is_high});
  }
  
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
    
    // 建立旧索引到新索引的映射
    std::vector<u32> old_to_new_mapping(current_pieces.size(), UINT32_MAX);
    
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
          
          // 记录映射关系
          old_to_new_mapping[i] = left_idx;
        } else {
          // 切割失败，保持原样
          size_t new_idx = new_pieces.size();
          new_pieces.push_back(piece);
          new_adjacency.resize(new_pieces.size());
          new_adjacency[new_idx] = current_adjacency[i];
          old_to_new_mapping[i] = new_idx;
        }
      } else {
        // 不需要切割
        size_t new_idx = new_pieces.size();
        new_pieces.push_back(piece);
        new_adjacency.resize(new_pieces.size());
        new_adjacency[new_idx] = current_adjacency[i];
        old_to_new_mapping[i] = new_idx;
      }
    }
    
    // 现在安全地建立邻接关系
    for(size_t i = 0; i < current_pieces.size(); ++i) {
      if(old_to_new_mapping[i] == UINT32_MAX) continue;
      
      const auto& piece = current_pieces[i];
      
      if(piece.bb.minx * 2 < line.coord && line.coord < piece.bb.maxx * 2) {
        // 这个片段被切割了
        auto cut_result = cut_polygon_with_adjacency(piece, false, line.coord, line.is_high);
        
        if(cut_result.pieces.size() == 2) {
          size_t left_idx = old_to_new_mapping[i];
          size_t right_idx = left_idx + 1;
          
          // 继承原片段的邻接关系
          for(u32 old_neighbor : current_adjacency[i]) {
            if(old_to_new_mapping[old_neighbor] != UINT32_MAX) {
              u32 new_neighbor = old_to_new_mapping[old_neighbor];
              new_adjacency[left_idx].push_back(new_neighbor);
              new_adjacency[new_neighbor].push_back(left_idx);
              new_adjacency[right_idx].push_back(new_neighbor);
              new_adjacency[new_neighbor].push_back(right_idx);
            }
          }
        }
      }
    }
    
    current_pieces = std::move(new_pieces);
    current_adjacency = std::move(new_adjacency);
  }
  
  // 再横切
  for(const auto& line : plan.horizontal_lines) {
    std::vector<Poly> new_pieces;
    std::vector<std::vector<u32>> new_adjacency;
    
    // 建立旧索引到新索引的映射
    std::vector<u32> old_to_new_mapping(current_pieces.size(), UINT32_MAX);
    
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
          
          // 记录映射关系
          old_to_new_mapping[i] = lower_idx;
        } else {
          // 切割失败，保持原样
          size_t new_idx = new_pieces.size();
          new_pieces.push_back(piece);
          new_adjacency.resize(new_pieces.size());
          new_adjacency[new_idx] = current_adjacency[i];
          old_to_new_mapping[i] = new_idx;
        }
      } else {
        // 不需要切割
        size_t new_idx = new_pieces.size();
        new_pieces.push_back(piece);
        new_adjacency.resize(new_pieces.size());
        new_adjacency[new_idx] = current_adjacency[i];
        old_to_new_mapping[i] = new_idx;
      }
    }
    
    // 现在安全地建立邻接关系
    for(size_t i = 0; i < current_pieces.size(); ++i) {
      if(old_to_new_mapping[i] == UINT32_MAX) continue;
      
      const auto& piece = current_pieces[i];
      
      if(piece.bb.miny * 2 < line.coord && line.coord < piece.bb.maxy * 2) {
        // 这个片段被切割了
        auto cut_result = cut_polygon_with_adjacency(piece, true, line.coord, line.is_high);
        
        if(cut_result.pieces.size() == 2) {
          size_t lower_idx = old_to_new_mapping[i];
          size_t upper_idx = lower_idx + 1;
          
          // 继承原片段的邻接关系
          for(u32 old_neighbor : current_adjacency[i]) {
            if(old_to_new_mapping[old_neighbor] != UINT32_MAX) {
              u32 new_neighbor = old_to_new_mapping[old_neighbor];
              new_adjacency[lower_idx].push_back(new_neighbor);
              new_adjacency[new_neighbor].push_back(lower_idx);
              new_adjacency[upper_idx].push_back(new_neighbor);
              new_adjacency[new_neighbor].push_back(upper_idx);
            }
          }
        }
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
  
  // 找到与u_gid对应的片段
  for(size_t i = 0; i < slices.pieces.size(); ++i) {
    if(slices.pieces[i].gid == u_gid) {
      u64 slice_gid = (u64(aa_lid) << 32) | (aa_pid << 16) | i;
      if(vis.find(slice_gid) == vis.end()) {
        q.push(slice_gid);
        vis.insert(slice_gid);
      }
    }
  }
}

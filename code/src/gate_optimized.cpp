#include "gate.hpp"
#include "geom.hpp"

#include <algorithm>
#include <unordered_set>
#include <vector>
#include <iostream>

namespace {

// ============================================================================
// 贯穿类型枚举
// ============================================================================
enum class PenetrationType {
  NONE,   // 无贯穿
  HORIZ,  // 横向贯穿：Gate横跨AA宽度，生成水平切割线
  VERT    // 纵向贯穿：Gate横跨AA高度，生成垂直切割线
};

// ============================================================================
// 切割计划结构
// ============================================================================
struct CutPlan {
  bool valid = false;           // 是否有效
  bool horizontal = false;      // 切割线方向：true=水平，false=垂直
  i32 coord = 0;               // 切割线坐标（整数单位）
  i32 span_min = 0;            // 切割线的跨度起点
  i32 span_max = 0;            // 切割线的跨度终点
  bool is_high = false;        // 电位属性
};

// ============================================================================
// BBox工具函数
// ============================================================================
BBox bbox_intersection(const BBox& a, const BBox& b) {
  BBox r;
  r.minx = std::max(a.minx, b.minx);
  r.maxx = std::min(a.maxx, b.maxx);
  r.miny = std::max(a.miny, b.miny);
  r.maxy = std::min(a.maxy, b.maxy);
  if(r.minx > r.maxx || r.miny > r.maxy) {
    r.minx = r.maxx = r.miny = r.maxy = 0;
  }
  return r;
}

// ============================================================================
// 核心判据：BBox中线跨越 + 边界容差
// ============================================================================
PenetrationType check_penetration(const Poly& gate, const Poly& aa) {
  // 步骤1：快速排除无重叠的情况
  if(!bbox_overlap(gate.bb, aa.bb)) {
    return PenetrationType::NONE;
  }
  
  // 步骤2：计算AA中心线
  i32 midx = (aa.bb.minx + aa.bb.maxx) / 2;
  i32 midy = (aa.bb.miny + aa.bb.maxy) / 2;
  
  // 步骤3：设置容差（以整数坐标为单位）
  const i32 tol = 1;  // 半格单位容差
  
  // 步骤4：判断贯穿方向
  bool crosses_vert_mid = (gate.bb.minx < midx && gate.bb.maxx > midx);
  bool crosses_horiz_mid = (gate.bb.miny < midy && gate.bb.maxy > midy);
  
  // 优先判断纵向贯穿（Gate横跨AA高度，生成垂直切割线）
  if(crosses_vert_mid) {
    // 确保Gate触及或接近AA的上下边界
    if(gate.bb.miny <= aa.bb.miny + tol && gate.bb.maxy >= aa.bb.maxy - tol) {
      return PenetrationType::VERT;
    }
  }
  
  // 判断横向贯穿（Gate横跨AA宽度，生成水平切割线）
  if(crosses_horiz_mid) {
    // 确保Gate触及或接近AA的左右边界
    if(gate.bb.minx <= aa.bb.minx + tol && gate.bb.maxx >= aa.bb.maxx - tol) {
      return PenetrationType::HORIZ;
    }
  }
  
  return PenetrationType::NONE;
}

// ============================================================================
// 切割计划生成函数
// ============================================================================
CutPlan compute_cut_plan(const Poly& gate, const Poly& aa, bool is_high) {
  CutPlan plan;
  plan.valid = false;
  
  PenetrationType type = check_penetration(gate, aa);
  
  if(type == PenetrationType::VERT) {
    // 纵向贯穿：生成垂直切割线
    plan.valid = true;
    plan.horizontal = false;
    
    // 计算交集范围
    BBox intersection = bbox_intersection(gate.bb, aa.bb);
    
    // 取交集区域的中心x坐标作为竖切线位置
    plan.coord = (intersection.minx + intersection.maxx) / 2;
    
    // 跨度为交集的y范围
    plan.span_min = intersection.miny;
    plan.span_max = intersection.maxy;
    plan.is_high = is_high;
    
  } else if(type == PenetrationType::HORIZ) {
    // 横向贯穿：生成水平切割线
    plan.valid = true;
    plan.horizontal = true;
    
    // 计算交集范围
    BBox intersection = bbox_intersection(gate.bb, aa.bb);
    
    // 取交集区域的中心y坐标作为横切线位置
    plan.coord = (intersection.miny + intersection.maxy) / 2;
    
    // 跨度为交集的x范围
    plan.span_min = intersection.minx;
    plan.span_max = intersection.maxx;
    plan.is_high = is_high;
  }
  
  return plan;
}

// ============================================================================
// 切割线合并辅助函数
// ============================================================================
struct PlannedLine {
  bool horizontal = false;
  i32 coord2 = 0;  // 2x坐标空间
  bool is_high = false;
  i32 span_min2 = 0;
  i32 span_max2 = 0;
};

bool intervals_overlap2(i32 a0, i32 a1, i32 b0, i32 b1) {
  return std::max(a0, b0) < std::min(a1, b1);
}

void merge_planned_lines(std::vector<PlannedLine>& lines) {
  // 排序：先按方向，再按坐标，再按跨度
  std::sort(lines.begin(), lines.end(), [](const PlannedLine& a, const PlannedLine& b) {
    if(a.horizontal != b.horizontal) return a.horizontal < b.horizontal;
    if(a.coord2 != b.coord2) return a.coord2 < b.coord2;
    if(a.span_min2 != b.span_min2) return a.span_min2 < b.span_min2;
    return a.span_max2 < b.span_max2;
  });

  // 合并重叠或相邻的线段
  std::vector<PlannedLine> merged;
  for(const auto& line : lines) {
    if(!merged.empty() &&
       merged.back().horizontal == line.horizontal &&
       merged.back().coord2 == line.coord2 &&
       intervals_overlap2(merged.back().span_min2, merged.back().span_max2,
                          line.span_min2, line.span_max2)) {
      // 合并：扩展跨度范围，电位取或
      merged.back().span_min2 = std::min(merged.back().span_min2, line.span_min2);
      merged.back().span_max2 = std::max(merged.back().span_max2, line.span_max2);
      merged.back().is_high = merged.back().is_high || line.is_high;
    } else {
      merged.push_back(line);
    }
  }
  lines.swap(merged);
}

// ============================================================================
// 多边形裁剪辅助函数（用于Lazy-Cut切割）
// ============================================================================
struct Pt2 {
  i32 x2 = 0;
  i32 y2 = 0;
};

std::vector<Pt2> to_pt2(const std::vector<Pt>& v) {
  std::vector<Pt2> r;
  r.reserve(v.size());
  for(const auto& p : v) {
    r.push_back({p.x * 2, p.y * 2});
  }
  return r;
}

Pt2 interpolate(const Pt2& a, const Pt2& b, i32 coord2, bool horizontal) {
  if(horizontal) {
    if(a.y2 == b.y2) return a;
    long long num = (long long)(coord2 - a.y2) * (b.x2 - a.x2);
    long long den = (long long)(b.y2 - a.y2);
    Pt2 r;
    r.y2 = coord2;
    if(den == 0) {
      r.x2 = a.x2;
    } else {
      r.x2 = a.x2 + (i32)(num / den);
    }
    return r;
  } else {
    if(a.x2 == b.x2) return a;
    long long num = (long long)(coord2 - a.x2) * (b.y2 - a.y2);
    long long den = (long long)(b.x2 - a.x2);
    Pt2 r;
    r.x2 = coord2;
    if(den == 0) {
      r.y2 = a.y2;
    } else {
      r.y2 = a.y2 + (i32)(num / den);
    }
    return r;
  }
}

std::vector<Pt2> clip_halfplane(const std::vector<Pt2>& subject,
                                bool horizontal,
                                i32 coord2,
                                bool keep_lower) {
  if(subject.empty()) return {};
  
  auto inside = [&](const Pt2& p) {
    if(horizontal) {
      return keep_lower ? (p.y2 <= coord2) : (p.y2 >= coord2);
    } else {
      return keep_lower ? (p.x2 <= coord2) : (p.x2 >= coord2);
    }
  };

  std::vector<Pt2> output;
  output.reserve(subject.size());

  Pt2 S = subject.back();
  bool S_inside = inside(S);
  for(const Pt2& E : subject) {
    bool E_inside = inside(E);
    if(E_inside) {
      if(!S_inside) {
        output.push_back(interpolate(S, E, coord2, horizontal));
      }
      output.push_back(E);
    } else if(S_inside) {
      output.push_back(interpolate(S, E, coord2, horizontal));
    }
    S = E;
    S_inside = E_inside;
  }

  // 去重
  std::vector<Pt2> dedup;
  dedup.reserve(output.size());
  for(const auto& p : output) {
    if(!dedup.empty() && dedup.back().x2 == p.x2 && dedup.back().y2 == p.y2) continue;
    dedup.push_back(p);
  }
  if(dedup.size() >= 2 && dedup.front().x2 == dedup.back().x2 && dedup.front().y2 == dedup.back().y2) {
    dedup.pop_back();
  }
  return dedup;
}

std::vector<Pt2> clip_polygon_axis(const std::vector<Pt2>& subject,
                                   bool horizontal,
                                   i32 coord2,
                                   bool keep_lower) {
  auto clipped = clip_halfplane(subject, horizontal, coord2, keep_lower);
  
  // 处理半格情况
  if(coord2 & 1) {
    for(auto& p : clipped) {
      if(horizontal) {
        if(p.y2 == coord2) p.y2 += keep_lower ? -1 : 1;
      } else {
        if(p.x2 == coord2) p.x2 += keep_lower ? -1 : 1;
      }
    }
  }
  
  if(clipped.size() < 3) return {};
  return clipped;
}

bool convert_vertices(const std::vector<Pt2>& v2, std::vector<Pt>& out) {
  out.clear();
  out.reserve(v2.size());
  for(const auto& p : v2) {
    if((p.x2 & 1) || (p.y2 & 1)) {
      return false;
    }
    out.push_back({p.x2 / 2, p.y2 / 2});
  }
  if(out.size() >= 2 && out.front().x == out.back().x && out.front().y == out.back().y) {
    out.pop_back();
  }
  return out.size() >= 3;
}

void rebuild_edges(Poly& poly) {
  poly.H.clear();
  poly.V.clear();
  const auto& v = poly.v;
  size_t n = v.size();
  for(size_t i = 0; i < n; ++i) {
    const Pt& a = v[i];
    const Pt& b = v[(i + 1) % n];
    if(a.x == b.x) {
      poly.V.push_back({a, b});
    } else if(a.y == b.y) {
      poly.H.push_back({a, b});
    }
  }
  poly.bb = bbox_of(poly.v);
}

} // namespace

// ============================================================================
// GateCtx 实现
// ============================================================================

void GateCtx::build_candidates(const Layout& L, int cell_hint) {
  if(poly_lid >= L.layers.size() || aa_lid >= L.layers.size()) return;

  const auto& poly_layer = L.layers[poly_lid];
  const auto& aa_layer = L.layers[aa_lid];

  poly_high.assign(poly_layer.polys.size(), 0);
  poly_to_aa_cands.assign(poly_layer.polys.size(), {});
  aa_to_poly_cands.assign(aa_layer.polys.size(), {});

  poly_grid = Grid::build(poly_layer.polys, cell_hint);
  aa_grid = Grid::build(aa_layer.polys, cell_hint);

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
  std::vector<PlannedLine> collected;

  // 遍历所有候选Gate多边形
  for(u32 poly_pid : aa_to_poly_cands[aa_pid]) {
    if(poly_pid >= poly_high.size()) continue;
    const auto& poly = L.layers[poly_lid].polys[poly_pid];
    
    // 使用新的判据计算切割计划
    bool is_high = (poly_high[poly_pid] != 0);
    CutPlan cut_plan = compute_cut_plan(poly, aa, is_high);
    
    if(cut_plan.valid) {
      // 转换为PlannedLine（2x坐标空间）
      PlannedLine line;
      line.horizontal = cut_plan.horizontal;
      line.coord2 = cut_plan.coord * 2;
      line.span_min2 = cut_plan.span_min * 2;
      line.span_max2 = cut_plan.span_max * 2;
      line.is_high = cut_plan.is_high;
      collected.push_back(line);
    }
  }

  // 合并重叠的切割线
  merge_planned_lines(collected);

  // 保存到plan
  plan.vertical_lines.clear();
  plan.horizontal_lines.clear();
  for(const auto& line : collected) {
    if(line.horizontal) {
      plan.horizontal_lines.push_back({line.coord2, line.is_high, line.span_min2, line.span_max2});
    } else {
      plan.vertical_lines.push_back({line.coord2, line.is_high, line.span_min2, line.span_max2});
    }
  }
  plan.computed = true;
}

void GateCtx::execute_cut(const Layout& L, u32 aa_pid) {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  
  auto& slices = aa_slices[aa_pid];
  if(slices.ready) return; // 已切割
  
  compute_aa_plan(L, aa_pid);
  
  const auto& plan = aa_plans[aa_pid];
  const auto& aa = L.layers[aa_lid].polys[aa_pid];

  // 如果没有切割线，返回原始AA
  if(plan.vertical_lines.empty() && plan.horizontal_lines.empty()) {
    slices.pieces.push_back(aa);
    slices.adjacency.resize(1);
    slices.ready = true;
    return;
  }

  // 执行切割：先处理垂直切割线
  std::vector<Poly> pieces = {aa};
  std::vector<std::vector<u32>> adjacency;
  std::vector<std::pair<u32,u32>> connections;
  
  u32 next_temp_id = 0;

  for(const auto& vline : plan.vertical_lines) {
    std::vector<Poly> new_pieces;
    std::vector<std::pair<u32,u32>> new_connections;
    
    for(size_t i = 0; i < pieces.size(); ++i) {
      const auto& piece = pieces[i];
      auto v2 = to_pt2(piece.v);
      
      // 左半部分
      auto left = clip_polygon_axis(v2, false, vline.coord, true);
      std::vector<Pt> left_v;
      bool has_left = convert_vertices(left, left_v);
      
      // 右半部分
      auto right = clip_polygon_axis(v2, false, vline.coord, false);
      std::vector<Pt> right_v;
      bool has_right = convert_vertices(right, right_v);
      
      if(has_left && has_right) {
        // 切割成功：左右两部分
        Poly left_poly;
        left_poly.v = left_v;
        left_poly.lid = aa.lid;  // 保持原AA的层ID
        rebuild_edges(left_poly);
        u32 left_id = new_pieces.size();
        new_pieces.push_back(left_poly);
        
        Poly right_poly;
        right_poly.v = right_v;
        right_poly.lid = aa.lid;  // 保持原AA的层ID
        rebuild_edges(right_poly);
        u32 right_id = new_pieces.size();
        new_pieces.push_back(right_poly);
        
        // 记录连接关系
        new_connections.push_back({left_id, right_id});
      } else if(has_left) {
        new_pieces.push_back(piece);
      } else if(has_right) {
        new_pieces.push_back(piece);
      }
    }
    
    pieces.swap(new_pieces);
    connections.insert(connections.end(), new_connections.begin(), new_connections.end());
  }

  // 处理水平切割线
  for(const auto& hline : plan.horizontal_lines) {
    std::vector<Poly> new_pieces;
    std::vector<std::pair<u32,u32>> new_connections;
    
    for(size_t i = 0; i < pieces.size(); ++i) {
      const auto& piece = pieces[i];
      auto v2 = to_pt2(piece.v);
      
      // 下半部分
      auto bottom = clip_polygon_axis(v2, true, hline.coord, true);
      std::vector<Pt> bottom_v;
      bool has_bottom = convert_vertices(bottom, bottom_v);
      
      // 上半部分
      auto top = clip_polygon_axis(v2, true, hline.coord, false);
      std::vector<Pt> top_v;
      bool has_top = convert_vertices(top, top_v);
      
      if(has_bottom && has_top) {
        // 切割成功：上下两部分
        Poly bottom_poly;
        bottom_poly.v = bottom_v;
        bottom_poly.lid = aa.lid;  // 保持原AA的层ID
        rebuild_edges(bottom_poly);
        u32 bottom_id = new_pieces.size();
        new_pieces.push_back(bottom_poly);
        
        Poly top_poly;
        top_poly.v = top_v;
        top_poly.lid = aa.lid;  // 保持原AA的层ID
        rebuild_edges(top_poly);
        u32 top_id = new_pieces.size();
        new_pieces.push_back(top_poly);
        
        // 记录连接关系
        new_connections.push_back({bottom_id, top_id});
      } else if(has_bottom) {
        new_pieces.push_back(piece);
      } else if(has_top) {
        new_pieces.push_back(piece);
      }
    }
    
    pieces.swap(new_pieces);
    connections.insert(connections.end(), new_connections.begin(), new_connections.end());
  }

  // 构建邻接表
  adjacency.resize(pieces.size());
  for(const auto& conn : connections) {
    adjacency[conn.first].push_back(conn.second);
    adjacency[conn.second].push_back(conn.first);
  }

  slices.pieces = pieces;
  slices.adjacency = adjacency;
  slices.ready = true;
}

void GateCtx::enqueue_entry_slices_from(u64 u_gid, u32 aa_pid, 
                                        std::queue<u64>& q, 
                                        std::unordered_set<u64>& vis,
                                        const Layout& L) const {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  
  // 确保AA已经切割
  auto slices_it = aa_slices.find(aa_pid);
  if(slices_it == aa_slices.end() || !slices_it->second.ready) {
    return;
  }
  
  const auto& slices = slices_it->second;
  
  // 获取源多边形（可能是普通多边形或切片）
  const Poly* source = nullptr;
  if(is_slice_gid(u_gid)) {
    // 如果u_gid本身就是一个切片
    source = slice_from_gid(u_gid);
  } else {
    // 如果是普通的多边形gid
    u32 u_lid = gid_lid(u_gid);
    u32 u_pid = gid_pid(u_gid);
    if(u_lid < L.layers.size() && u_pid < L.layers[u_lid].polys.size()) {
      source = &L.layers[u_lid].polys[u_pid];
    }
  }
  
  if(!source) return;
  
  // 找到与源多边形相交的AA片段（入口片段）
  for(size_t i = 0; i < slices.pieces.size(); ++i) {
    const auto& slice = slices.pieces[i];
    
    // 快速BBox检查
    if(!bbox_overlap(source->bb, slice.bb)) continue;
    
    // 精确的多边形相交检查
    if(!poly_intersect_manhattan(*source, slice)) continue;
    
    // 这是一个入口片段，加入队列
    u64 slice_gid = make_slice_gid(aa_pid, i);
    if(vis.find(slice_gid) == vis.end()) {
      vis.insert(slice_gid);
      q.push(slice_gid);
      
      // 沿着导通边在AA内部扩展（高电平切割线连接的片段）
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


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

// 计算曼哈顿多边形的面积（Shoelace公式）
int64_t compute_polygon_area(const std::vector<Pt>& v) {
  if(v.size() < 3) return 0;
  int64_t area2 = 0;  // 2倍面积
  for(size_t i = 0; i < v.size(); ++i) {
    const auto& p1 = v[i];
    const auto& p2 = v[(i + 1) % v.size()];
    area2 += (int64_t)p1.x * p2.y - (int64_t)p2.x * p1.y;
  }
  return std::abs(area2) / 2;
}

// ============================================================================
// 核心判据：简单宽松判定（只要有实质性相交就切割）
// 几何切割与电平无关：所有相交的AA都要切，高/低电平只影响切后导通性
// ============================================================================
PenetrationType check_penetration(const Poly& gate, const Poly& aa) {
  // 步骤1：bbox预筛
  if(!bbox_overlap(gate.bb, aa.bb)) {
    return PenetrationType::NONE;
  }
  
  // 步骤2：计算交集区域 I = Gate ∩ AA
  BBox I;
  I.minx = std::max(gate.bb.minx, aa.bb.minx);
  I.miny = std::max(gate.bb.miny, aa.bb.miny);
  I.maxx = std::min(gate.bb.maxx, aa.bb.maxx);
  I.maxy = std::min(gate.bb.maxy, aa.bb.maxy);
  
  if(I.minx >= I.maxx || I.miny >= I.maxy) {
    return PenetrationType::NONE;  // 无有效交集
  }
  
  // 步骤3：简单判定 - 只要有实质性相交就切割
  i32 I_width = I.maxx - I.minx;
  i32 I_height = I.maxy - I.miny;
  
  // 最小相交阈值（避免擦边误切）
  const i32 min_intersection = 2;
  
  if(I_width < min_intersection || I_height < min_intersection) {
    return PenetrationType::NONE;  // 交集太小，视为擦边
  }
  
  // 决定切割方向：根据Gate的形状选择
  i32 gate_width = gate.bb.maxx - gate.bb.minx;
  i32 gate_height = gate.bb.maxy - gate.bb.miny;
  
  // 如果Gate是竖长的（高>宽），用垂直切割
  // 如果Gate是横长的（宽>高），用水平切割
  if(gate_height > gate_width) {
    return PenetrationType::VERT;
  } else {
    return PenetrationType::HORIZ;
  }
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

bool intervals_overlap_or_adjacent(i32 a0, i32 a1, i32 b0, i32 b1) {
  // 检查两个区间是否重叠或相邻（可以合并）
  // 重叠: [a0,a1] 和 [b0,b1] 有公共部分
  // 相邻: a1 == b0 或 b1 == a0
  return std::max(a0, b0) <= std::min(a1, b1);
}

void merge_planned_lines(std::vector<PlannedLine>& lines) {
  if(lines.empty()) return;
  
  size_t original_count = lines.size();
  
  // 步骤1: 坐标对齐（snap to grid），减少微小偏差
  const i32 SNAP = 2;  // 2x坐标空间的对齐步长
  auto snap = [](i32 c2) -> i32 {
    return ((c2 + SNAP/2) / SNAP) * SNAP;
  };
  
  int snap_changes = 0;
  for(auto& ln : lines) {
    i32 old_coord = ln.coord2;
    ln.coord2 = snap(ln.coord2);
    if(old_coord != ln.coord2) snap_changes++;
  }
  
  // 步骤2: 排序：先按方向，再按坐标，再按跨度
  std::sort(lines.begin(), lines.end(), [](const PlannedLine& a, const PlannedLine& b) {
    if(a.horizontal != b.horizontal) return a.horizontal < b.horizontal;
    if(a.coord2 != b.coord2) return a.coord2 < b.coord2;
    if(a.span_min2 != b.span_min2) return a.span_min2 < b.span_min2;
    return a.span_max2 < b.span_max2;
  });

  // 步骤3: 合并重叠或相邻的线段（改进3：按电平分桶，严禁高低电平合并）
  const i32 COORD_TOL = 4;   // 允许坐标偏差（2x空间）
  const i32 SPAN_TOL = 2;    // 跨度相邻容差
  
  int merge_count = 0;
  std::vector<PlannedLine> merged;
  for(const auto& line : lines) {
    if(!merged.empty()) {
      auto& m = merged.back();
      bool same_dir = (m.horizontal == line.horizontal);
      bool same_coord = (std::abs(m.coord2 - line.coord2) <= COORD_TOL);
      bool same_polarity = (m.is_high == line.is_high);  // 改进3：电平必须一致
      
      // 检查跨度是否重叠或相邻（放宽容差）
      i32 gap = std::max(m.span_min2, line.span_min2) - std::min(m.span_max2, line.span_max2);
      bool span_meet = (gap <= SPAN_TOL);
      
      if(same_dir && same_coord && same_polarity && span_meet) {
        // 合并：扩展跨度范围（电平已保证一致，无需OR）
        m.span_min2 = std::min(m.span_min2, line.span_min2);
        m.span_max2 = std::max(m.span_max2, line.span_max2);
        // m.is_high 保持不变（已经一致）
        merge_count++;
        continue;
      }
    }
    merged.push_back(line);
  }
  
  lines.swap(merged);
  
  // 调试输出（仅在有合并时）
  if(merge_count > 0 || snap_changes > 0) {
    std::cerr << "[merge_planned_lines] original=" << original_count 
              << " snap_changes=" << snap_changes
              << " merged=" << merge_count
              << " final=" << lines.size() << "\n";
  }
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
  // 确保顶点逆时针顺序
  ensure_ccw(poly.v);
  
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
// 辅助函数：检查两个切片是否相邻并且应该导通
// ============================================================================

// 检查两个片段是否共享一条边界，并且该边界对应的切割线是高电平
static bool are_adjacent_slices(const Poly& slice1, const Poly& slice2, const GateCtx::AAPlan& plan) {
  // 快速排除：BBox不相邻
  if(!bbox_overlap(slice1.bb, slice2.bb)) {
    return false;
  }
  
  // 检查是否共享边界
  // 两个片段相邻意味着它们有一条共享的切割线
  
  // 方法：检查是否有顶点在对方的边界上
  // 这里采用简化版本：检查BBox是否在某个维度上紧邻
  
  const i32 tol = 2;  // 2x坐标空间的容差
  
  // 检查垂直切割线（左右相邻）
  bool share_vertical_edge = false;
  i32 shared_x = 0;
  
  if(std::abs(slice1.bb.maxx - slice2.bb.minx) <= tol) {
    // slice1在左，slice2在右
    share_vertical_edge = true;
    shared_x = (slice1.bb.maxx + slice2.bb.minx) / 2;
  } else if(std::abs(slice2.bb.maxx - slice1.bb.minx) <= tol) {
    // slice2在左，slice1在右
    share_vertical_edge = true;
    shared_x = (slice2.bb.maxx + slice1.bb.minx) / 2;
  }
  
  if(share_vertical_edge) {
    // 检查这条垂直切割线是否是高电平
    for(const auto& vline : plan.vertical_lines) {
      if(std::abs(vline.coord - shared_x) <= tol) {
        // 找到对应的切割线，检查是否高电平
        if(vline.is_high) {
          // 还需要检查y方向是否有重叠
          i32 y_overlap_min = std::max(slice1.bb.miny, slice2.bb.miny);
          i32 y_overlap_max = std::min(slice1.bb.maxy, slice2.bb.maxy);
          if(y_overlap_max > y_overlap_min) {
            return true;  // 相邻且导通
          }
        }
        break;
      }
    }
  }
  
  // 检查水平切割线（上下相邻）
  bool share_horizontal_edge = false;
  i32 shared_y = 0;
  
  if(std::abs(slice1.bb.maxy - slice2.bb.miny) <= tol) {
    // slice1在下，slice2在上
    share_horizontal_edge = true;
    shared_y = (slice1.bb.maxy + slice2.bb.miny) / 2;
  } else if(std::abs(slice2.bb.maxy - slice1.bb.miny) <= tol) {
    // slice2在下，slice1在上
    share_horizontal_edge = true;
    shared_y = (slice2.bb.maxy + slice1.bb.miny) / 2;
  }
  
  if(share_horizontal_edge) {
    // 检查这条水平切割线是否是高电平
    for(const auto& hline : plan.horizontal_lines) {
      if(std::abs(hline.coord - shared_y) <= tol) {
        if(hline.is_high) {
          // 检查x方向是否有重叠
          i32 x_overlap_min = std::max(slice1.bb.minx, slice2.bb.minx);
          i32 x_overlap_max = std::min(slice1.bb.maxx, slice2.bb.maxx);
          if(x_overlap_max > x_overlap_min) {
            return true;
          }
        }
        break;
      }
    }
  }
  
  return false;
}

// ============================================================================
// 切割线驱动的邻接边建立（优化B）
// ============================================================================

// 对于垂直切割线：找到所有与该线相邻的片段，分为左右两组，并连接
static void connect_by_vertical_line(
    const GateCtx::CutLine& line,
    const std::vector<Poly>& pieces,
    std::vector<std::vector<u32>>& adjacency) {
  
  if(!line.is_high) return;  // 只有高电平线才建立导通
  
  static int call_count = 0;
  call_count++;
  bool debug = (call_count <= 5);  // 只调试前5次调用
  
  const i32 tol = 10;  // 放宽容差从4到10，允许更大的对齐误差
  // 重要：line的坐标是2x空间，需要转换为1x空间来匹配pieces的bbox
  i32 coord = line.coord / 2;
  i32 y_min = line.span_min2 / 2;
  i32 y_max = line.span_max2 / 2;
  
  if(debug) {
    std::cerr << "[connect_by_vertical_line #" << call_count << "] coord=" << coord 
              << " (2x=" << line.coord << ") y_range=[" << y_min << "," << y_max << "] pieces=" << pieces.size() << "\n";
  }
  
  std::vector<u32> left_pieces, right_pieces;
  
  for(u32 i = 0; i < pieces.size(); ++i) {
    const auto& bb = pieces[i].bb;
    
    // 检查片段是否与切割线在x方向对齐
    bool touches_line = (bb.minx <= coord + tol && bb.maxx >= coord - tol);
    if(!touches_line) continue;
    
    // 检查y方向是否有重叠（接触也算，放宽判定）
    i32 y_overlap_min = std::max(bb.miny, y_min);
    i32 y_overlap_max = std::min(bb.maxy, y_max);
    if(y_overlap_max < y_overlap_min - tol) continue;  // 改为允许接触
    
    // 判断片段在切割线的左侧还是右侧
    if(bb.maxx <= coord + tol) {
      left_pieces.push_back(i);
      if(debug) std::cerr << "  piece[" << i << "] on LEFT, bbox=[" << bb.minx << "," << bb.miny << "," << bb.maxx << "," << bb.maxy << "]\n";
    } else if(bb.minx >= coord - tol) {
      right_pieces.push_back(i);
      if(debug) std::cerr << "  piece[" << i << "] on RIGHT, bbox=[" << bb.minx << "," << bb.miny << "," << bb.maxx << "," << bb.maxy << "]\n";
    }
  }
  
  if(debug) {
    std::cerr << "  Found " << left_pieces.size() << " left pieces, " << right_pieces.size() << " right pieces\n";
  }
  
  // 连接左右两侧y区间重叠的片段对
  int connections_made = 0;
  for(u32 li : left_pieces) {
    for(u32 ri : right_pieces) {
      const auto& lbb = pieces[li].bb;
      const auto& rbb = pieces[ri].bb;
      
      // 检查y区间是否重叠（接触也算，放宽判定）
      i32 y_overlap_min = std::max(lbb.miny, rbb.miny);
      i32 y_overlap_max = std::min(lbb.maxy, rbb.maxy);
      
      if(y_overlap_max >= y_overlap_min - tol) {  // 改为允许接触
        // 建立双向边
        adjacency[li].push_back(ri);
        adjacency[ri].push_back(li);
        connections_made++;
        if(debug) std::cerr << "  Connected piece[" << li << "] <-> piece[" << ri << "]\n";
      }
    }
  }
  
  if(debug) {
    std::cerr << "  Total connections made: " << connections_made << "\n";
  }
}

// 对于水平切割线：找到所有与该线相邻的片段，分为上下两组，并连接
static void connect_by_horizontal_line(
    const GateCtx::CutLine& line,
    const std::vector<Poly>& pieces,
    std::vector<std::vector<u32>>& adjacency) {
  
  if(!line.is_high) return;  // 只有高电平线才建立导通
  
  const i32 tol = 10;  // 放宽容差从4到10，允许更大的对齐误差
  // 重要：line的坐标是2x空间，需要转换为1x空间来匹配pieces的bbox
  i32 coord = line.coord / 2;
  i32 x_min = line.span_min2 / 2;
  i32 x_max = line.span_max2 / 2;
  
  std::vector<u32> bottom_pieces, top_pieces;
  
  for(u32 i = 0; i < pieces.size(); ++i) {
    const auto& bb = pieces[i].bb;
    
    // 检查片段是否与切割线在y方向对齐
    bool touches_line = (bb.miny <= coord + tol && bb.maxy >= coord - tol);
    if(!touches_line) continue;
    
    // 检查x方向是否有重叠（接触也算，放宽判定）
    i32 x_overlap_min = std::max(bb.minx, x_min);
    i32 x_overlap_max = std::min(bb.maxx, x_max);
    if(x_overlap_max < x_overlap_min - tol) continue;  // 改为允许接触
    
    // 判断片段在切割线的下方还是上方
    if(bb.maxy <= coord + tol) {
      bottom_pieces.push_back(i);
    } else if(bb.miny >= coord - tol) {
      top_pieces.push_back(i);
    }
  }
  
  // 连接上下两侧x区间重叠的片段对
  for(u32 bi : bottom_pieces) {
    for(u32 ti : top_pieces) {
      const auto& bbb = pieces[bi].bb;
      const auto& tbb = pieces[ti].bb;
      
      // 检查x区间是否重叠（接触也算，放宽判定）
      i32 x_overlap_min = std::max(bbb.minx, tbb.minx);
      i32 x_overlap_max = std::min(bbb.maxx, tbb.maxx);
      
      if(x_overlap_max >= x_overlap_min - tol) {  // 改为允许接触
        // 建立双向边
        adjacency[bi].push_back(ti);
        adjacency[ti].push_back(bi);
      }
    }
  }
}

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

  // 改进4：候选对齐升级到"确实相交"，降低"高电平但不贯穿"噪声
  int bbox_only = 0, real_intersect = 0;
  for(u32 poly_pid = 0; poly_pid < poly_layer.polys.size(); ++poly_pid) {
    const auto& poly = poly_layer.polys[poly_pid];
    std::vector<u32> candidates;
    aa_grid.query(poly.bb, candidates);
    for(u32 aa_pid : candidates) {
      const auto& aa = aa_layer.polys[aa_pid];
      if(!bbox_overlap(poly.bb, aa.bb)) continue;
      bbox_only++;
      
      // 升级：使用几何相交判断代替纯bbox
      if(poly_intersect_manhattan(poly, aa)) {
        poly_to_aa_cands[poly_pid].push_back(aa_pid);
        aa_to_poly_cands[aa_pid].push_back(poly_pid);
        real_intersect++;
      }
    }
  }
  
  std::cerr << "[build_candidates] bbox_overlap=" << bbox_only 
            << ", real_intersect=" << real_intersect 
            << " (" << (bbox_only > 0 ? real_intersect*100/bbox_only : 0) << "%)\n";
}

void GateCtx::compute_aa_plan(const Layout& L, u32 aa_pid) {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  auto& plan = aa_plans[aa_pid];
  if(plan.computed) return;

  const auto& aa = L.layers[aa_lid].polys[aa_pid];
  std::vector<PlannedLine> collected;

  // 全局统计：分析未切AA的原因
  static int total_aa_processed = 0;
  static int aa_no_candidates = 0;  // 没有任何候选Poly
  static int aa_has_cands_but_no_penetrate = 0;  // 有候选但都不贯穿
  static int aa_will_be_cut = 0;  // 将被切割
  
  static int debug_count = 0;
  static int high_found = 0, low_found = 0;
  static int checked_high_but_not_penetrate = 0;
  
  total_aa_processed++;
  
  if(aa_pid == 0) {
    std::cerr << "[compute_aa_plan aa_pid=0] aa_to_poly_cands.size()=" << aa_to_poly_cands[aa_pid].size() << "\n";
    // 检查有多少poly_high被标记了
    int marked_count = 0;
    for(size_t i = 0; i < poly_high.size(); ++i) {
      if(poly_high[i] != 0) marked_count++;
    }
    std::cerr << "[compute_aa_plan aa_pid=0] poly_high marked count=" << marked_count << " / " << poly_high.size() << "\n";
  }
  
  // 检查是否有候选
  if(aa_to_poly_cands[aa_pid].empty()) {
    aa_no_candidates++;
  }
  
  bool has_any_penetration = false;
  
  for(u32 poly_pid : aa_to_poly_cands[aa_pid]) {
    // 优化C：确保只处理PO层的多边形
    if(poly_pid >= L.layers[poly_lid].polys.size()) continue;
    if(poly_pid >= poly_high.size()) continue;
    const auto& poly = L.layers[poly_lid].polys[poly_pid];
    // 二次确认：poly必须在poly_lid层
    if(poly.lid != poly_lid) {
      std::cerr << "[WARNING] poly_pid=" << poly_pid << " has lid=" << poly.lid 
                << " but expected poly_lid=" << poly_lid << "\n";
      continue;
    }
    
    // 使用新的判据计算切割计划
    bool is_high = (poly_high[poly_pid] != 0);
    CutPlan cut_plan = compute_cut_plan(poly, aa, is_high);
    
    if(cut_plan.valid) {
      has_any_penetration = true;
      
      if(is_high) high_found++;
      else low_found++;
      
      if(debug_count < 5) {
        std::cerr << "[compute_aa_plan] poly_pid=" << poly_pid 
                  << " is_high=" << is_high 
                  << " (poly_high[" << poly_pid << "]=" << (int)poly_high[poly_pid] << ")\n";
        debug_count++;
      }
      
      // 转换为PlannedLine（2x坐标空间）
      PlannedLine line;
      line.horizontal = cut_plan.horizontal;
      line.coord2 = cut_plan.coord * 2;
      line.span_min2 = cut_plan.span_min * 2;
      line.span_max2 = cut_plan.span_max * 2;
      line.is_high = cut_plan.is_high;
      collected.push_back(line);
    } else if(is_high) {
      // 高电平但不贯穿
      checked_high_but_not_penetrate++;
    }
  }
  
  // 统计未切原因
  if(!has_any_penetration && !aa_to_poly_cands[aa_pid].empty()) {
    aa_has_cands_but_no_penetrate++;
  }
  
  if(has_any_penetration) {
    aa_will_be_cut++;
  }
  
  if(aa_pid == 0) {
    std::cerr << "[compute_aa_plan] Total penetrations: high=" << high_found << " low=" << low_found << "\n";
    std::cerr << "[compute_aa_plan] High polys that don't penetrate: " << checked_high_but_not_penetrate << "\n";
  }
  
  // 打印最终统计（在最后一个AA处理完后）
  if(aa_pid == L.layers[aa_lid].polys.size() - 1) {
    std::cerr << "\n========================================\n";
    std::cerr << "  未切AA诊断报告\n";
    std::cerr << "========================================\n";
    std::cerr << "总AA数: " << total_aa_processed << "\n";
    std::cerr << "将被切割: " << aa_will_be_cut << " (" << (aa_will_be_cut*100/total_aa_processed) << "%)\n";
    std::cerr << "不切割: " << (total_aa_processed - aa_will_be_cut) << " (" << ((total_aa_processed-aa_will_be_cut)*100/total_aa_processed) << "%)\n";
    std::cerr << "\n未切原因分析:\n";
    std::cerr << "  1. 没有候选Poly: " << aa_no_candidates << "\n";
    std::cerr << "  2. 有候选但都不贯穿: " << aa_has_cands_but_no_penetrate << "\n";
    std::cerr << "  总计未切: " << (aa_no_candidates + aa_has_cands_but_no_penetrate) << "\n";
    std::cerr << "\n贯穿统计:\n";
    std::cerr << "  高电平贯穿: " << high_found << "\n";
    std::cerr << "  低电平贯穿: " << low_found << "\n";
    std::cerr << "  高电平但不贯穿: " << checked_high_but_not_penetrate << "\n";
    std::cerr << "========================================\n\n";
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
  
  // 调试特定AA
  bool debug_this_aa = (aa_pid == 26584);
  if(debug_this_aa) {
    std::cerr << "\n=== 调试AA pid=" << aa_pid << " ===\n";
    std::cerr << "AA bbox: [" << aa.bb.minx << "," << aa.bb.miny << "] 到 [" 
              << aa.bb.maxx << "," << aa.bb.maxy << "]\n";
    std::cerr << "切割线数量: vert=" << plan.vertical_lines.size() 
              << ", horiz=" << plan.horizontal_lines.size() << "\n";
    for(size_t i = 0; i < plan.vertical_lines.size(); i++) {
      const auto& vl = plan.vertical_lines[i];
      std::cerr << "  竖线#" << i << ": x=" << vl.coord/2 << " (2x=" << vl.coord 
                << "), y∈[" << vl.span_min2/2 << "," << vl.span_max2/2 << "], "
                << (vl.is_high ? "高电平" : "低电平") << "\n";
    }
    for(size_t i = 0; i < plan.horizontal_lines.size(); i++) {
      const auto& hl = plan.horizontal_lines[i];
      std::cerr << "  横线#" << i << ": y=" << hl.coord/2 << " (2x=" << hl.coord 
                << "), x∈[" << hl.span_min2/2 << "," << hl.span_max2/2 << "], "
                << (hl.is_high ? "高电平" : "低电平") << "\n";
    }
  }

  // 如果没有切割线，返回原始AA
  if(plan.vertical_lines.empty() && plan.horizontal_lines.empty()) {
    slices.pieces.push_back(aa);
    slices.adjacency.resize(1);
    slices.ready = true;
    return;
  }

  // 执行切割：先处理垂直切割线
  // 使用映射表追踪piece的演变关系
  std::vector<Poly> pieces = {aa};
  std::vector<std::vector<u32>> adjacency;
  
  // old_to_new[i] = 旧piece[i]对应的新piece索引列表（可能被切成多个）
  std::vector<std::vector<u32>> old_to_new;
  old_to_new.push_back({0}); // 初始：原始AA映射到piece[0]

  for(const auto& vline : plan.vertical_lines) {
    std::vector<Poly> new_pieces;
    std::vector<std::vector<u32>> new_old_to_new;
    
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
        left_poly.lid = aa.lid;
        rebuild_edges(left_poly);
        u32 left_id = new_pieces.size();
        new_pieces.push_back(left_poly);
        
        Poly right_poly;
        right_poly.v = right_v;
        right_poly.lid = aa.lid;
        rebuild_edges(right_poly);
        u32 right_id = new_pieces.size();
        new_pieces.push_back(right_poly);
        
        // 记录映射：旧piece[i] -> {left_id, right_id}
        new_old_to_new.push_back({left_id, right_id});
      } else if(has_left) {
        u32 new_id = new_pieces.size();
        new_pieces.push_back(piece);
        new_old_to_new.push_back({new_id});
      } else if(has_right) {
        u32 new_id = new_pieces.size();
        new_pieces.push_back(piece);
        new_old_to_new.push_back({new_id});
      } else {
        // piece完全被裁掉（不应该发生）
        new_old_to_new.push_back({});
      }
    }
    
    pieces.swap(new_pieces);
    old_to_new.swap(new_old_to_new);
  }

  // 处理水平切割线
  for(const auto& hline : plan.horizontal_lines) {
    std::vector<Poly> new_pieces;
    std::vector<std::vector<u32>> new_old_to_new;
    
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
        bottom_poly.lid = aa.lid;
        rebuild_edges(bottom_poly);
        u32 bottom_id = new_pieces.size();
        new_pieces.push_back(bottom_poly);
        
        Poly top_poly;
        top_poly.v = top_v;
        top_poly.lid = aa.lid;
        rebuild_edges(top_poly);
        u32 top_id = new_pieces.size();
        new_pieces.push_back(top_poly);
        
        new_old_to_new.push_back({bottom_id, top_id});
      } else if(has_bottom) {
        u32 new_id = new_pieces.size();
        new_pieces.push_back(piece);
        new_old_to_new.push_back({new_id});
      } else if(has_top) {
        u32 new_id = new_pieces.size();
        new_pieces.push_back(piece);
        new_old_to_new.push_back({new_id});
      } else {
        new_old_to_new.push_back({});
      }
    }
    
    pieces.swap(new_pieces);
    old_to_new.swap(new_old_to_new);
  }

  // 重新建立邻接关系：使用切割线驱动的方法（优化B）
  adjacency.resize(pieces.size());
  
  static int total_edges = 0;
  static bool first_aa = true;
  
  int local_edges = 0;
  
  // 遍历所有垂直切割线，建立左右相邻片段的导通边
  for(const auto& vline : plan.vertical_lines) {
    size_t before = total_edges;
    connect_by_vertical_line(vline, pieces, adjacency);
    size_t added = 0;
    for(const auto& adj_list : adjacency) {
      added += adj_list.size();
    }
    local_edges += (added - before);
  }
  
  // 遍历所有水平切割线，建立上下相邻片段的导通边
  for(const auto& hline : plan.horizontal_lines) {
    size_t before = 0;
    for(const auto& adj_list : adjacency) {
      before += adj_list.size();
    }
    connect_by_horizontal_line(hline, pieces, adjacency);
    size_t after = 0;
    for(const auto& adj_list : adjacency) {
      after += adj_list.size();
    }
    local_edges += (after - before);
  }
  
  // 统计实际邻接边数量（每条边被计数两次）
  size_t total_adj_entries = 0;
  for(const auto& adj_list : adjacency) {
    total_adj_entries += adj_list.size();
  }
  int actual_edges = total_adj_entries / 2;
  
  // 统计切线分布（前10个AA）
  static int aa_count = 0;
  aa_count++;
  
  if(aa_count <= 10 && pieces.size() > 1) {
    // 统计高电平切割线数量
    int high_vlines = 0, low_vlines = 0;
    int high_hlines = 0, low_hlines = 0;
    for(const auto& vl : plan.vertical_lines) {
      if(vl.is_high) high_vlines++;
      else low_vlines++;
    }
    for(const auto& hl : plan.horizontal_lines) {
      if(hl.is_high) high_hlines++;
      else low_hlines++;
    }
    
    std::cerr << "[AA #" << aa_count << "] aa_pid=" << aa_pid 
              << " pieces=" << pieces.size() 
              << " vlines=" << plan.vertical_lines.size() << " (high=" << high_vlines << ", low=" << low_vlines << ")"
              << " hlines=" << plan.horizontal_lines.size() << " (high=" << high_hlines << ", low=" << low_hlines << ")"
              << " adjacency_edges=" << actual_edges << "\n";
    
    // 打印所有切片的bbox
    for(size_t i = 0; i < pieces.size(); ++i) {
      const auto& p = pieces[i];
      std::cerr << "  piece[" << i << "] bbox=[" << p.bb.minx << "," << p.bb.miny << "," << p.bb.maxx << "," << p.bb.maxy << "]\n";
    }
    
    // 打印所有高电平切线（coord和span都是2x坐标）
    for(const auto& vl : plan.vertical_lines) {
      if(vl.is_high) {
        std::cerr << "  high_vline: coord=" << vl.coord << " (1x=" << vl.coord/2 << ") span2=[" << vl.span_min2 << "," << vl.span_max2 << "]\n";
      }
    }
    first_aa = false;
  }
  
  total_edges += actual_edges;

  slices.pieces = pieces;
  slices.adjacency = adjacency;
  slices.ready = true;
  
  if(debug_this_aa) {
    std::cerr << "\n切割完成！生成 " << pieces.size() << " 个片段：\n";
    for(size_t i = 0; i < pieces.size(); i++) {
      const auto& p = pieces[i];
      std::cerr << "  片段#" << i << ": bbox=[" << p.bb.minx << "," << p.bb.miny 
                << "] 到 [" << p.bb.maxx << "," << p.bb.maxy << "], "
                << "中心x=" << (p.bb.minx+p.bb.maxx)/2 << ", "
                << "邻居数=" << adjacency[i].size() << "\n";
      if(!adjacency[i].empty()) {
        std::cerr << "    邻居: ";
        for(u32 nb : adjacency[i]) {
          std::cerr << "#" << nb << " ";
        }
        std::cerr << "\n";
      }
    }
  }
  
  // 全局统计高低电平切线
  static int total_high_vlines = 0;
  static int total_low_vlines = 0;
  static int total_high_hlines = 0;
  static int total_low_hlines = 0;
  static bool stats_printed = false;
  
  for(const auto& vl : plan.vertical_lines) {
    if(vl.is_high) total_high_vlines++;
    else total_low_vlines++;
  }
  for(const auto& hl : plan.horizontal_lines) {
    if(hl.is_high) total_high_hlines++;
    else total_low_hlines++;
  }
  
  // 在最后打印统计
  static int aa_cut_count = 0;
  aa_cut_count++;
  if(aa_cut_count >= 26400 && !stats_printed) {
    stats_printed = true;
    int total_high = total_high_vlines + total_high_hlines;
    int total_low = total_low_vlines + total_low_hlines;
    int total_all = total_high + total_low;
    std::cerr << "\n=== Global Cutting Line Statistics ===\n";
    std::cerr << "High-potential vertical lines: " << total_high_vlines << "\n";
    std::cerr << "Low-potential vertical lines: " << total_low_vlines << "\n";
    std::cerr << "High-potential horizontal lines: " << total_high_hlines << "\n";
    std::cerr << "Low-potential horizontal lines: " << total_low_hlines << "\n";
    std::cerr << "Total high-potential lines: " << total_high << " (" 
              << (total_all > 0 ? (100.0 * total_high / total_all) : 0) << "%)\n";
    std::cerr << "Total low-potential lines: " << total_low << " (" 
              << (total_all > 0 ? (100.0 * total_low / total_all) : 0) << "%)\n";
    std::cerr << "======================================\n\n";
  }
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


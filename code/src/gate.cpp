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

// 更严格的贯穿判据：检查poly是否真正"穿过"AA
bool is_full_penetration(const Poly& poly, const Poly& aa) {
  BBox intersection = bbox_intersection(poly.bb, aa.bb);
  
  // 检查是否真正贯穿
  // 横向贯穿：交集必须接触AA的上下边界，且poly在AA的左右两侧都有延伸
  bool horizontal_penetration = (intersection.miny <= aa.bb.miny && 
                                intersection.maxy >= aa.bb.maxy);
  
  // 纵向贯穿：交集必须接触AA的左右边界，且poly在AA的上下两侧都有延伸
  bool vertical_penetration = (intersection.minx <= aa.bb.minx && 
                              intersection.maxx >= aa.bb.maxx);
  
  // 进一步验证：检查poly是否真的"穿过"AA
  if(horizontal_penetration) {
    // 检查poly是否在AA的左右两侧都有延伸
    bool extends_left = (poly.bb.minx < aa.bb.minx);
    bool extends_right = (poly.bb.maxx > aa.bb.maxx);
    if(!extends_left || !extends_right) return false;
  }
  
  if(vertical_penetration) {
    // 检查poly是否在AA的上下两侧都有延伸
    bool extends_bottom = (poly.bb.miny < aa.bb.miny);
    bool extends_top = (poly.bb.maxy > aa.bb.maxy);
    if(!extends_bottom || !extends_top) return false;
  }
  
  return horizontal_penetration || vertical_penetration;
}

// 沿轴对齐线切割多边形（纯整数算法，支持半格切线）
// is_horizontal: true表示横切（y=coord2/2），false表示竖切（x=coord2/2）
// coord2: 2×坐标（支持半格）
// keep_lower: true保留下侧/左侧，false保留上侧/右侧
std::vector<Pt> clip_axis_aligned_halfgrid(const std::vector<Pt>& vertices, bool is_horizontal, i32 coord2, bool keep_lower) {
  std::vector<Pt> result;
  int n = vertices.size();
  if(n == 0) return result;
  
  for(int i = 0; i < n; ++i) {
    const Pt& current = vertices[i];
    const Pt& next = vertices[(i + 1) % n];
    
    // 获取切割维度的坐标（2×）
    i32 curr_cut2 = is_horizontal ? (current.y * 2) : (current.x * 2);
    i32 next_cut2 = is_horizontal ? (next.y * 2) : (next.x * 2);
    
    // 获取另一维度的坐标
    i32 curr_other = is_horizontal ? current.x : current.y;
    i32 next_other = is_horizontal ? next.x : next.y;
    
    // 判断点是否在保留侧（用2×坐标判断）
    bool curr_inside = keep_lower ? (curr_cut2 < coord2) : (curr_cut2 >= coord2);
    bool next_inside = keep_lower ? (next_cut2 < coord2) : (next_cut2 >= coord2);
    
    if(curr_inside) {
      result.push_back(current);
    }
    
    // 边与切线相交
    if(curr_inside != next_inside) {
      // 边 current->next 穿越切线
      if(curr_cut2 == next_cut2) {
        // 边平行于切线
        continue;
      }
      
      // 计算交点
      // 如果coord2是偶数（整数切线），交点是整数
      // 如果coord2是奇数（半格切线），交点在水平/垂直边上
      
      if(is_horizontal) {
        // 横切：y方向的边（垂直边）
        if(current.x == next.x) {
          // 垂直边，交点x坐标就是current.x
          // y坐标：coord2是奇数时不需要输出（落在边内部）
          if(coord2 % 2 == 0) {
            result.push_back({current.x, coord2 / 2});
          }
          // 半格情况：S-H会在两个端点处理，这里不输出
        } else {
          // 水平边，根据coord2是否在边的y范围内
          // 实际上水平边不应该穿越横切线（因为curr_inside != next_inside意味着y不同）
        }
      } else {
        // 竖切：x方向的边（水平边）
        if(current.y == next.y) {
          // 水平边，交点y坐标就是current.y
          if(coord2 % 2 == 0) {
            result.push_back({coord2 / 2, current.y});
          }
        } else {
          // 垂直边
        }
      }
    }
  }
  
  return result;
}

// 计算多边形在水平线y=const处的x区间集合（内部覆盖的部分）
std::vector<std::pair<i32, i32>> xIntervalsOnHorizontal(const Poly& poly, i32 y) {
  std::vector<i32> x_crossings;
  
  for(size_t i = 0; i < poly.v.size(); ++i) {
    const Pt& p1 = poly.v[i];
    const Pt& p2 = poly.v[(i + 1) % poly.v.size()];
    
    i32 y1 = p1.y, y2 = p2.y;
    i32 x1 = p1.x, x2 = p2.x;
    
    // 边与y线相交
    if((y1 <= y && y < y2) || (y2 <= y && y < y1)) {
      // 垂直边
      if(x1 == x2) {
        x_crossings.push_back(x1);
      } else {
        // 斜边（曼哈顿多边形不应有，但保险起见）
        long long x_at_y = x1 + (long long)(y - y1) * (x2 - x1) / (y2 - y1);
        x_crossings.push_back((i32)x_at_y);
      }
    } else if(y1 == y && y2 == y) {
      // 水平边恰好在y线上：两个端点都算
      x_crossings.push_back(x1);
      x_crossings.push_back(x2);
    }
  }
  
  // 升序并配对
  std::sort(x_crossings.begin(), x_crossings.end());
  std::vector<std::pair<i32, i32>> intervals;
  for(size_t i = 0; i + 1 < x_crossings.size(); i += 2) {
    intervals.push_back({x_crossings[i], x_crossings[i+1]});
  }
  
  return intervals;
}

// 计算多边形在垂直线x=const处的y区间集合
std::vector<std::pair<i32, i32>> yIntervalsOnVertical(const Poly& poly, i32 x) {
  std::vector<i32> y_crossings;
  
  for(size_t i = 0; i < poly.v.size(); ++i) {
    const Pt& p1 = poly.v[i];
    const Pt& p2 = poly.v[(i + 1) % poly.v.size()];
    
    i32 x1 = p1.x, x2 = p2.x;
    i32 y1 = p1.y, y2 = p2.y;
    
    if((x1 <= x && x < x2) || (x2 <= x && x < x1)) {
      if(y1 == y2) {
        y_crossings.push_back(y1);
      } else {
        long long y_at_x = y1 + (long long)(x - x1) * (y2 - y1) / (x2 - x1);
        y_crossings.push_back((i32)y_at_x);
      }
    } else if(x1 == x && x2 == x) {
      y_crossings.push_back(y1);
      y_crossings.push_back(y2);
    }
  }
  
  std::sort(y_crossings.begin(), y_crossings.end());
  std::vector<std::pair<i32, i32>> intervals;
  for(size_t i = 0; i + 1 < y_crossings.size(); i += 2) {
    intervals.push_back({y_crossings[i], y_crossings[i+1]});
  }
  
  return intervals;
}

} // namespace

void GateCtx::build_candidates(const Layout& L, int cell_hint){
  if(poly_lid >= L.layers.size() || aa_lid >= L.layers.size()) return;
  
  const auto& poly_layer = L.layers[poly_lid];
  const auto& aa_layer = L.layers[aa_lid];
  
  poly_grid = Grid::build(poly_layer.polys, cell_hint);
  aa_grid = Grid::build(aa_layer.polys, cell_hint);
  
  aa_to_poly_cands.resize(aa_layer.polys.size());
  for(u32 aa_pid=0; aa_pid<aa_layer.polys.size(); ++aa_pid){
    const auto& AA = aa_layer.polys[aa_pid];
    std::vector<u32> cands;
    poly_grid.query(AA.bb, cands);
    std::sort(cands.begin(), cands.end());
    cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
    aa_to_poly_cands[aa_pid] = std::move(cands);
  }
}

void GateCtx::compute_aa_plan(const Layout& L, u32 aa_pid){
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  if(aa_plans[aa_pid].computed) return;
  
  const Poly& AA = L.layers[aa_lid].polys[aa_pid];
  AAPlan& plan = aa_plans[aa_pid];
  
  // 使用2×坐标系，支持半格切线
  std::unordered_set<i32> vertical_lines2_high, vertical_lines2_low;
  std::unordered_set<i32> horizontal_lines2_high, horizontal_lines2_low;
  
  // 直接用AA.bb查询Poly层，确保全覆盖
  std::vector<u32> cands;
  poly_grid.query(AA.bb, cands);
  std::sort(cands.begin(), cands.end());
  cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
  
  for(u32 poly_pid : cands) {
    if(poly_pid >= L.layers[poly_lid].polys.size()) continue;
    
    const Poly& PP = L.layers[poly_lid].polys[poly_pid];
    
    // 粗筛：bbox重叠
    if(!bbox_overlap(AA.bb, PP.bb)) continue;
    
    // 精判：曼哈顿几何相交
    if(!poly_intersect_manhattan(AA, PP)) continue;
    
    // 更严格的贯穿判据：检查是否真正贯穿
    if(!is_full_penetration(PP, AA)) continue;
    
    bool is_high = (poly_pid < poly_high.size() && poly_high[poly_pid] != 0);
    
    // 竖向贯穿：区间交法
    auto topSegs = xIntervalsOnHorizontal(PP, AA.bb.maxy);
    auto bottomSegs = xIntervalsOnHorizontal(PP, AA.bb.miny);
    
    for(const auto& [lT, rT] : topSegs) {
      for(const auto& [lB, rB] : bottomSegs) {
        i32 L = std::max(lT, lB);
        i32 R = std::min(rT, rB);
        if(L < R) {
          // 通道存在，生成半格切线坐标
          i32 cut2 = L + R;  // 2×中线坐标
          
          // 过滤AA边界外的
          if(cut2 <= AA.bb.minx * 2 || cut2 >= AA.bb.maxx * 2) continue;
          
          // 高/低分类（同坐标高优先）
          if(is_high) {
            vertical_lines2_high.insert(cut2);
            vertical_lines2_low.erase(cut2);
          } else if(vertical_lines2_high.find(cut2) == vertical_lines2_high.end()) {
            vertical_lines2_low.insert(cut2);
          }
        }
      }
    }
    
    // 横向贯穿：区间交法
    auto leftSegs = yIntervalsOnVertical(PP, AA.bb.minx);
    auto rightSegs = yIntervalsOnVertical(PP, AA.bb.maxx);
    
    for(const auto& [lL, rL] : leftSegs) {
      for(const auto& [lR, rR] : rightSegs) {
        i32 L = std::max(lL, lR);
        i32 R = std::min(rL, rR);
        if(L < R) {
          i32 cut2 = L + R;
          
          if(cut2 <= AA.bb.miny * 2 || cut2 >= AA.bb.maxy * 2) continue;
          
          if(is_high) {
            horizontal_lines2_high.insert(cut2);
            horizontal_lines2_low.erase(cut2);
          } else if(horizontal_lines2_high.find(cut2) == horizontal_lines2_high.end()) {
            horizontal_lines2_low.insert(cut2);
          }
        }
      }
    }
  }
  
  // 转换为有序列表（仍用2×坐标）
  for(i32 cut2 : vertical_lines2_high) {
    plan.vertical_lines.push_back({cut2, true});
  }
  for(i32 cut2 : vertical_lines2_low) {
    plan.vertical_lines.push_back({cut2, false});
  }
  for(i32 cut2 : horizontal_lines2_high) {
    plan.horizontal_lines.push_back({cut2, true});
  }
  for(i32 cut2 : horizontal_lines2_low) {
    plan.horizontal_lines.push_back({cut2, false});
  }
  
  // 排序
  std::sort(plan.vertical_lines.begin(), plan.vertical_lines.end(), 
            [](const CutLine& a, const CutLine& b) { return a.coord < b.coord; });
  std::sort(plan.horizontal_lines.begin(), plan.horizontal_lines.end(),
            [](const CutLine& a, const CutLine& b) { return a.coord < b.coord; });
  
  plan.computed = true;
}

void GateCtx::execute_cut(const Layout& L, u32 aa_pid){
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  if(aa_slices[aa_pid].ready) return;
  
  // 确保切割计划已计算
  compute_aa_plan(L, aa_pid);
  
  const Poly& AA = L.layers[aa_lid].polys[aa_pid];
  const AAPlan& plan = aa_plans[aa_pid];
  
  // 如果没有切割线，保持原样
  if(plan.vertical_lines.empty() && plan.horizontal_lines.empty()) {
    AASlices& slices = aa_slices[aa_pid];
    slices.pieces.push_back(AA);
    slices.pieces[0].pid = 0;
    slices.pieces[0].gid = make_slice_gid(aa_pid, 0);
    slices.adjacency.resize(1);
    slices.ready = true;
    return;
  }
  
  // 提取切割线坐标
  std::vector<i32> v_coords, h_coords;
  for(const auto& line : plan.vertical_lines) {
    v_coords.push_back(line.coord);
  }
  for(const auto& line : plan.horizontal_lines) {
    h_coords.push_back(line.coord);
  }
  
  // 网格化切割与导通边重建
  
  // 1. 构造有序坐标（包含AA边界和所有切线，转换为单×坐标）
  std::set<i32> X_set, Y_set;
  X_set.insert(AA.bb.minx);
  X_set.insert(AA.bb.maxx);
  
  for(const auto& line : plan.vertical_lines) {
    // 2×坐标转换：如果是偶数，是整数切线；如果是奇数，是半格切线
    if(line.coord % 2 == 0) {
      X_set.insert(line.coord / 2);
    } else {
      // 半格切线：在(coord/2)和(coord/2+1)之间，添加两个边界
      X_set.insert(line.coord / 2);
      X_set.insert(line.coord / 2 + 1);
    }
  }
  
  for(const auto& line : plan.horizontal_lines) {
    if(line.coord % 2 == 0) {
      Y_set.insert(line.coord / 2);
    } else {
      Y_set.insert(line.coord / 2);
      Y_set.insert(line.coord / 2 + 1);
    }
  }
  
  std::vector<i32> X(X_set.begin(), X_set.end());
  std::vector<i32> Y(Y_set.begin(), Y_set.end());
  
  // 2. 记录哪些2×坐标是高电平切线
  std::unordered_set<i32> high_vertical2, high_horizontal2;
  for(const auto& line : plan.vertical_lines) {
    if(line.is_high) high_vertical2.insert(line.coord);
  }
  for(const auto& line : plan.horizontal_lines) {
    if(line.is_high) high_horizontal2.insert(line.coord);
  }
  
  // 3. 执行切割（先竖后横）
  std::vector<Poly> pieces;
  pieces.push_back(AA);
  
  // 先竖切（使用2×坐标，支持半格切线）
  for(const auto& line : plan.vertical_lines) {
    i32 x_cut2 = line.coord;  // 2×坐标
    std::vector<Poly> new_pieces;
    
    for(const auto& piece : pieces) {
      // 判断是否需要切割（用2×坐标）
      if(piece.bb.minx * 2 < x_cut2 && x_cut2 < piece.bb.maxx * 2) {
        // 切割这个片段
        auto left = clip_axis_aligned_halfgrid(piece.v, false, x_cut2, true);
        auto right = clip_axis_aligned_halfgrid(piece.v, false, x_cut2, false);
        
        if(left.size() >= 3) {
          Poly p_left;
          p_left.v = std::move(left);
          p_left.bb = bbox_of(p_left.v);
          p_left.lid = piece.lid;
          new_pieces.push_back(std::move(p_left));
        }
        if(right.size() >= 3) {
          Poly p_right;
          p_right.v = std::move(right);
          p_right.bb = bbox_of(p_right.v);
          p_right.lid = piece.lid;
          new_pieces.push_back(std::move(p_right));
        }
      } else {
        new_pieces.push_back(piece);
      }
    }
    
    pieces = std::move(new_pieces);
  }
  
  // 再横切（使用2×坐标，支持半格切线）
  for(const auto& line : plan.horizontal_lines) {
    i32 y_cut2 = line.coord;  // 2×坐标
    std::vector<Poly> new_pieces;
    
    for(const auto& piece : pieces) {
      if(piece.bb.miny * 2 < y_cut2 && y_cut2 < piece.bb.maxy * 2) {
        // 切割这个片段
        auto lower = clip_axis_aligned_halfgrid(piece.v, true, y_cut2, true);
        auto upper = clip_axis_aligned_halfgrid(piece.v, true, y_cut2, false);
        
        if(lower.size() >= 3) {
          Poly p_lower;
          p_lower.v = std::move(lower);
          p_lower.bb = bbox_of(p_lower.v);
          p_lower.lid = piece.lid;
          new_pieces.push_back(std::move(p_lower));
        }
        if(upper.size() >= 3) {
          Poly p_upper;
          p_upper.v = std::move(upper);
          p_upper.bb = bbox_of(p_upper.v);
          p_upper.lid = piece.lid;
          new_pieces.push_back(std::move(p_upper));
        }
      } else {
        new_pieces.push_back(piece);
      }
    }
    
    pieces = std::move(new_pieces);
  }
  
  // 4. 建立网格单元映射：piece ⇔ cell(i,j)
  std::unordered_map<u64, u32> cell2piece;  // (i<<32)|j -> piece_idx
  int mapped_count = 0;
  int unmapped_count = 0;
  
  for(size_t piece_idx = 0; piece_idx < pieces.size(); ++piece_idx) {
    const Poly& piece = pieces[piece_idx];
    
    // 查找piece在网格中的位置
    auto ix0_it = std::lower_bound(X.begin(), X.end(), piece.bb.minx);
    auto ix1_it = std::lower_bound(X.begin(), X.end(), piece.bb.maxx);
    auto iy0_it = std::lower_bound(Y.begin(), Y.end(), piece.bb.miny);
    auto iy1_it = std::lower_bound(Y.begin(), Y.end(), piece.bb.maxy);
    
    if(ix0_it == X.end() || ix1_it == X.end() || 
       iy0_it == Y.end() || iy1_it == Y.end()) {
      unmapped_count++;
      continue;
    }
    
    int ix0 = ix0_it - X.begin();
    int ix1 = ix1_it - X.begin();
    int iy0 = iy0_it - Y.begin();
    int iy1 = iy1_it - Y.begin();
    
    // 验证：片段bbox应该等于[X[i], X[i+1]] × [Y[j], Y[j+1]]
    if(X[ix0] == piece.bb.minx && X[ix1] == piece.bb.maxx && ix1 == ix0 + 1 &&
       Y[iy0] == piece.bb.miny && Y[iy1] == piece.bb.maxy && iy1 == iy0 + 1) {
      // 完美匹配网格单元
      u64 cell_key = ((u64)ix0 << 32) | iy0;
      cell2piece[cell_key] = piece_idx;
      mapped_count++;
    } else {
      unmapped_count++;
    }
  }
  
  if(aa_pid == 0) {  // 只打印第一个AA的调试信息
    std::cerr << "Debug AA[" << aa_pid << "]: " << pieces.size() << " pieces, "
              << "X.size=" << X.size() << ", Y.size=" << Y.size() << ", "
              << "mapped=" << mapped_count << ", unmapped=" << unmapped_count << "\n";
    std::cerr << "  Vertical lines: " << plan.vertical_lines.size() 
              << " (high: " << high_vertical2.size() << ")\n";
    std::cerr << "  Horizontal lines: " << plan.horizontal_lines.size()
              << " (high: " << high_horizontal2.size() << ")\n";
    std::cerr << "  X coords: ";
    for(auto x : X) std::cerr << x << " ";
    std::cerr << "\n  Y coords: ";
    for(auto y : Y) std::cerr << y << " ";
    std::cerr << "\n";
  }
  
  // 5. 网格化重建导通边
  std::vector<std::vector<u32>> adjacency(pieces.size());
  int vertical_edges = 0, horizontal_edges = 0;
  
  // 5a. 垂直高电平切线导通：检查X[k]处是否有高电平切线
  for(int k = 1; k < (int)X.size() - 1; ++k) {
    // 检查X[k]对应的2×坐标是否是高电平切线
    i32 x2 = X[k] * 2;
    bool is_high = (high_vertical2.find(x2) != high_vertical2.end());
    
    // 也检查半格切线：x2-1和x2+1
    if(!is_high && high_vertical2.find(x2-1) != high_vertical2.end()) is_high = true;
    if(!is_high && high_vertical2.find(x2+1) != high_vertical2.end()) is_high = true;
    
    if(!is_high) continue;
    
    // 枚举所有纵向间隔
    for(int j = 0; j < (int)Y.size() - 1; ++j) {
      int i_left = k - 1;
      int i_right = k;
      
      u64 cell_left = ((u64)i_left << 32) | j;
      u64 cell_right = ((u64)i_right << 32) | j;
      
      auto it_left = cell2piece.find(cell_left);
      auto it_right = cell2piece.find(cell_right);
      
      if(it_left != cell2piece.end() && it_right != cell2piece.end()) {
        u32 a = it_left->second;
        u32 b = it_right->second;
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
        vertical_edges++;
      }
    }
  }
  
  // 5b. 水平高电平切线导通
  for(int l = 1; l < (int)Y.size() - 1; ++l) {
    i32 y2 = Y[l] * 2;
    bool is_high = (high_horizontal2.find(y2) != high_horizontal2.end());
    
    if(!is_high && high_horizontal2.find(y2-1) != high_horizontal2.end()) is_high = true;
    if(!is_high && high_horizontal2.find(y2+1) != high_horizontal2.end()) is_high = true;
    
    if(!is_high) continue;
    
    // 枚举所有横向间隔
    for(int i = 0; i < (int)X.size() - 1; ++i) {
      int j_down = l - 1;
      int j_up = l;
      
      u64 cell_down = ((u64)i << 32) | j_down;
      u64 cell_up = ((u64)i << 32) | j_up;
      
      auto it_down = cell2piece.find(cell_down);
      auto it_up = cell2piece.find(cell_up);
      
      if(it_down != cell2piece.end() && it_up != cell2piece.end()) {
        u32 a = it_down->second;
        u32 b = it_up->second;
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
        horizontal_edges++;
      }
    }
  }
  
  if(aa_pid == 0) {  // 打印第一个AA的导通边统计
    std::cerr << "  Conductive edges: vertical=" << vertical_edges 
              << ", horizontal=" << horizontal_edges << "\n";
  }
  
  
  // 6. 去重导通边
  for(auto& neighbors : adjacency) {
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
  }
  
  // 7. 赋值gid和pid
  for(size_t i = 0; i < pieces.size(); ++i) {
    pieces[i].pid = i;
    pieces[i].gid = make_slice_gid(aa_pid, i);
  }
  
  AASlices& slices = aa_slices[aa_pid];
  slices.pieces = std::move(pieces);
  slices.adjacency = std::move(adjacency);
  slices.ready = true;
}

// 计算bbox重叠面积
static long long overlap_area(const BBox& a, const BBox& b) {
  i32 x_overlap = std::max(0, std::min(a.maxx, b.maxx) - std::max(a.minx, b.minx));
  i32 y_overlap = std::max(0, std::min(a.maxy, b.maxy) - std::max(a.miny, b.miny));
  return (long long)x_overlap * y_overlap;
}

void GateCtx::enqueue_entry_slices_from(u64 u_gid, u32 aa_pid,
                                        std::queue<u64>& q,
                                        std::unordered_set<u64>& vis,
                                        const Layout& L) const {
  // 确保AA已切割
  auto it_slices = aa_slices.find(aa_pid);
  if(it_slices == aa_slices.end() || !it_slices->second.ready) return;
  
  const auto& S = it_slices->second;
  const auto& pieces = S.pieces;
  
  // 获取U的多边形
  const Poly* U_poly = nullptr;
  if(is_slice_gid(u_gid)) {
    U_poly = slice_from_gid(u_gid);
  } else {
    u32 u_lid = gid_lid(u_gid);
    u32 u_pid = gid_pid(u_gid);
    if(u_lid < L.layers.size() && u_pid < L.layers[u_lid].polys.size()) {
      U_poly = &L.layers[u_lid].polys[u_pid];
    }
  }
  
  if(!U_poly) return;
  
  // 找到所有入口片（与U几何相交的片段）
  std::vector<int> entries;
  for(size_t i = 0; i < pieces.size(); ++i) {
    if(!bbox_overlap(pieces[i].bb, U_poly->bb)) continue;
    if(poly_intersect_manhattan(pieces[i], *U_poly)) {
      entries.push_back(i);
    }
  }
  
  // 兜底：如果没有精判相交，选bbox重叠面积最大的
  if(entries.empty()) {
    int best = -1;
    long long bestA = -1;
    for(size_t i = 0; i < pieces.size(); ++i) {
      long long A = overlap_area(pieces[i].bb, U_poly->bb);
      if(A > bestA) {
        bestA = A;
        best = i;
      }
    }
    if(best >= 0) entries.push_back(best);
  }
  
  // 在AA内部沿高切线邻接泛洪
  std::queue<int> iq;
  std::vector<char> seen(pieces.size(), 0);
  
  for(int sidx : entries) {
    u64 gid = make_slice_gid(aa_pid, sidx);
    if(vis.insert(gid).second) {
      q.push(gid);
      iq.push(sidx);
      seen[sidx] = 1;
    }
  }
  
  // 沿导通边（只走高切线建立的边）扩展
  while(!iq.empty()) {
    int u = iq.front(); iq.pop();
    if(u >= (int)S.adjacency.size()) continue;
    
    for(u32 v : S.adjacency[u]) {
      if(v >= pieces.size()) continue;
      u64 gid = make_slice_gid(aa_pid, v);
      if(!seen[v] && vis.insert(gid).second) {
        q.push(gid);
        iq.push(v);
        seen[v] = 1;
      }
    }
  }
}

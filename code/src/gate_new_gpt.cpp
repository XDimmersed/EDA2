#include "gate.hpp"
#include "geom.hpp"

#include <algorithm>
#include <unordered_set>
#include <vector>

namespace {

struct Pt2 {
  i32 x2 = 0;
  i32 y2 = 0;
};

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

struct PlannedLine {
  bool horizontal = false;
  i32 coord2 = 0;
  bool is_high = false;
  i32 span_min = 0;
  i32 span_max = 0;
};

bool intervals_overlap(i32 a0, i32 a1, i32 b0, i32 b1) {
  return std::max(a0, b0) < std::min(a1, b1);
}

void merge_planned_lines(std::vector<PlannedLine>& lines) {
  std::sort(lines.begin(), lines.end(), [](const PlannedLine& a, const PlannedLine& b) {
    if(a.horizontal != b.horizontal) return b.horizontal < a.horizontal;
    if(a.coord2 != b.coord2) return a.coord2 < b.coord2;
    if(a.span_min != b.span_min) return a.span_min < b.span_min;
    return a.span_max < b.span_max;
  });

  std::vector<PlannedLine> merged;
  for(const auto& line : lines) {
    if(!merged.empty() &&
       merged.back().horizontal == line.horizontal &&
       merged.back().coord2 == line.coord2 &&
       intervals_overlap(merged.back().span_min, merged.back().span_max,
                         line.span_min, line.span_max)) {
      merged.back().span_min = std::min(merged.back().span_min, line.span_min);
      merged.back().span_max = std::max(merged.back().span_max, line.span_max);
      merged.back().is_high = merged.back().is_high || line.is_high;
    } else {
      merged.push_back(line);
    }
  }
  lines.swap(merged);
}

} // namespace

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

  for(u32 poly_pid : aa_to_poly_cands[aa_pid]) {
    if(poly_pid >= poly_high.size()) continue;
    const auto& poly = L.layers[poly_lid].polys[poly_pid];
    if(!bbox_overlap(poly.bb, aa.bb)) continue;

    auto subject = to_pt2(poly.v);
    if(subject.size() < 3) continue;

    std::vector<Pt2> tmp = subject;
    tmp = clip_halfplane(tmp, false, aa.bb.minx * 2, false);
    tmp = clip_halfplane(tmp, false, aa.bb.maxx * 2, true);
    tmp = clip_halfplane(tmp, true, aa.bb.miny * 2, false);
    tmp = clip_halfplane(tmp, true, aa.bb.maxy * 2, true);
    if(tmp.size() < 3) continue;

    BBox inter_bb = {tmp[0].x2 / 2, tmp[0].y2 / 2, tmp[0].x2 / 2, tmp[0].y2 / 2};
    for(const auto& p : tmp) {
      inter_bb.minx = std::min(inter_bb.minx, p.x2 / 2);
      inter_bb.maxx = std::max(inter_bb.maxx, p.x2 / 2);
      inter_bb.miny = std::min(inter_bb.miny, p.y2 / 2);
      inter_bb.maxy = std::max(inter_bb.maxy, p.y2 / 2);
    }

    bool touches_left = (inter_bb.minx <= aa.bb.minx);
    bool touches_right = (inter_bb.maxx >= aa.bb.maxx);
    bool touches_bottom = (inter_bb.miny <= aa.bb.miny);
    bool touches_top = (inter_bb.maxy >= aa.bb.maxy);

    bool is_high = (poly_high[poly_pid] != 0);

    if(touches_left && touches_right && inter_bb.maxy > inter_bb.miny) {
      PlannedLine line;
      line.horizontal = true;
      line.coord2 = inter_bb.miny + inter_bb.maxy;
      line.span_min = inter_bb.minx;
      line.span_max = inter_bb.maxx;
      line.is_high = is_high;
      collected.push_back(line);
    }
    if(touches_bottom && touches_top && inter_bb.maxx > inter_bb.minx) {
      PlannedLine line;
      line.horizontal = false;
      line.coord2 = inter_bb.minx + inter_bb.maxx;
      line.span_min = inter_bb.miny;
      line.span_max = inter_bb.maxy;
      line.is_high = is_high;
      collected.push_back(line);
    }
  }

  merge_planned_lines(collected);

  plan.vertical_lines.clear();
  plan.horizontal_lines.clear();
  for(const auto& line : collected) {
    if(line.horizontal) {
      plan.horizontal_lines.push_back({line.coord2, line.is_high});
    } else {
      plan.vertical_lines.push_back({line.coord2, line.is_high});
    }
  }
  plan.computed = true;
}

static bool convert_poly_from_clip(const std::vector<Pt2>& v2, const Poly& base, Poly& out) {
  std::vector<Pt> verts;
  if(!convert_vertices(v2, verts)) return false;
  ensure_ccw(verts);
  out = base;
  out.v = std::move(verts);
  rebuild_edges(out);
  return true;
}

static std::vector<PlannedLine> gather_planned_lines(const GateCtx::AAPlan& plan,
                                                     const Poly& aa) {
  std::vector<PlannedLine> lines;
  for(const auto& line : plan.vertical_lines) {
    PlannedLine pl;
    pl.horizontal = false;
    pl.coord2 = line.coord;
    pl.is_high = line.is_high;
    pl.span_min = aa.bb.miny;
    pl.span_max = aa.bb.maxy;
    lines.push_back(pl);
  }
  for(const auto& line : plan.horizontal_lines) {
    PlannedLine pl;
    pl.horizontal = true;
    pl.coord2 = line.coord;
    pl.is_high = line.is_high;
    pl.span_min = aa.bb.minx;
    pl.span_max = aa.bb.maxx;
    lines.push_back(pl);
  }
  merge_planned_lines(lines);
  return lines;
}

void GateCtx::execute_cut(const Layout& L, u32 aa_pid) {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  auto& slices = aa_slices[aa_pid];
  if(slices.ready) return;

  compute_aa_plan(L, aa_pid);
  const auto& plan = aa_plans[aa_pid];
  const auto& aa = L.layers[aa_lid].polys[aa_pid];

  auto lines = gather_planned_lines(plan, aa);
  if(lines.empty()) {
    slices.pieces = {aa};
    slices.adjacency.assign(1, {});
    slices.ready = true;
    return;
  }

  std::vector<Poly> pieces = {aa};

  for(const auto& line : lines) {
    std::vector<Poly> next;
    next.reserve(pieces.size() * 2);
    for(const auto& piece : pieces) {
      if(line.horizontal) {
        if(piece.bb.miny * 2 >= line.coord2 || piece.bb.maxy * 2 <= line.coord2) {
          next.push_back(piece);
          continue;
        }
      } else {
        if(piece.bb.minx * 2 >= line.coord2 || piece.bb.maxx * 2 <= line.coord2) {
          next.push_back(piece);
          continue;
        }
      }

      auto subject = to_pt2(piece.v);
      auto lower = clip_polygon_axis(subject, line.horizontal, line.coord2, true);
      auto upper = clip_polygon_axis(subject, line.horizontal, line.coord2, false);
      Poly lower_poly, upper_poly;
      bool ok_lower = convert_poly_from_clip(lower, aa, lower_poly);
      bool ok_upper = convert_poly_from_clip(upper, aa, upper_poly);

      if(ok_lower) next.push_back(std::move(lower_poly));
      if(ok_upper) next.push_back(std::move(upper_poly));
      if(!ok_lower && !ok_upper) {
        next.push_back(piece);
      }
    }
    pieces.swap(next);
  }

  slices.pieces = pieces;
  slices.adjacency.assign(pieces.size(), {});

  auto register_pairs = [&](const PlannedLine& line) {
    std::vector<u32> side_a;
    std::vector<u32> side_b;
    for(u32 idx = 0; idx < slices.pieces.size(); ++idx) {
      const auto& p = slices.pieces[idx];
      if(line.horizontal) {
        if(p.bb.maxy * 2 == line.coord2) {
          side_a.push_back(idx);
        } else if(p.bb.miny * 2 == line.coord2) {
          side_b.push_back(idx);
        }
      } else {
        if(p.bb.maxx * 2 == line.coord2) {
          side_a.push_back(idx);
        } else if(p.bb.minx * 2 == line.coord2) {
          side_b.push_back(idx);
        }
      }
    }

    if(!line.is_high) return;

    for(u32 a : side_a) {
      for(u32 b : side_b) {
        const auto& pa = slices.pieces[a];
        const auto& pb = slices.pieces[b];
        if(line.horizontal) {
          if(!intervals_overlap(pa.bb.minx, pa.bb.maxx, pb.bb.minx, pb.bb.maxx)) continue;
          if(!intervals_overlap(pa.bb.minx, pa.bb.maxx, line.span_min, line.span_max)) continue;
        } else {
          if(!intervals_overlap(pa.bb.miny, pa.bb.maxy, pb.bb.miny, pb.bb.maxy)) continue;
          if(!intervals_overlap(pa.bb.miny, pa.bb.maxy, line.span_min, line.span_max)) continue;
        }
        slices.adjacency[a].push_back(b);
        slices.adjacency[b].push_back(a);
      }
    }
  };

  for(const auto& line : lines) {
    register_pairs(line);
  }

  for(auto& nbrs : slices.adjacency) {
    std::sort(nbrs.begin(), nbrs.end());
    nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
  }

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

  const Poly* source = nullptr;
  if(is_slice_gid(u_gid)) {
    source = slice_from_gid(u_gid);
  } else {
    source = &L.layers[gid_lid(u_gid)].polys[gid_pid(u_gid)];
  }
  if(!source) return;

  for(u32 idx = 0; idx < slices.pieces.size(); ++idx) {
    const auto& piece = slices.pieces[idx];
    if(!bbox_overlap(piece.bb, source->bb)) continue;
    if(!poly_intersect_manhattan(piece, *source)) continue;
    u64 gid = make_slice_gid(aa_pid, idx);
    if(vis.insert(gid).second) {
      q.push(gid);
      for(u32 nb : slices.adjacency[idx]) {
        u64 nb_gid = make_slice_gid(aa_pid, nb);
        if(vis.insert(nb_gid).second) {
          q.push(nb_gid);
        }
      }
    }
  }
}


#include "gate.hpp"
#include "geom.hpp"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct Rect {
  i32 x0 = 0, x1 = 0, y0 = 0, y1 = 0;
  bool high = false;
};

struct Boundary {
  bool open = false;
  bool blocked = false;
  void mark(bool is_high) {
    if(is_high) open = true;
    else blocked = true;
  }
};

Poly make_rect_poly(const Poly& base, i32 x0, i32 y0, i32 x1, i32 y1) {
  Poly out;
  out.lid = base.lid;
  out.bb = {x0, y0, x1, y1};
  out.v = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  out.H = { { {x0, y0}, {x1, y0} }, { {x1, y1}, {x0, y1} } };
  out.V = { { {x1, y0}, {x1, y1} }, { {x0, y1}, {x0, y0} } };
  return out;
}

struct CellRect {
  i32 x0 = 0, x1 = 0, y0 = 0, y1 = 0;
};

struct PointKey {
  i32 x = 0, y = 0;
  bool operator==(const PointKey& o) const { return x == o.x && y == o.y; }
  bool operator<(const PointKey& o) const {
    if(y != o.y) return y < o.y;
    return x < o.x;
  }
};

struct PointKeyHash {
  size_t operator()(const PointKey& p) const {
    return (static_cast<size_t>(p.x) << 32) ^ static_cast<size_t>(p.y);
  }
};

struct EdgeKey {
  PointKey a;
  PointKey b;
  bool operator==(const EdgeKey& o) const { return a.x == o.a.x && a.y == o.a.y && b.x == o.b.x && b.y == o.b.y; }
};

struct EdgeKeyHash {
  size_t operator()(const EdgeKey& e) const {
    return PointKeyHash{}(e.a) ^ (PointKeyHash{}(e.b) << 1);
  }
};

void add_edge(std::unordered_map<EdgeKey, int, EdgeKeyHash>& edges, const PointKey& a, const PointKey& b) {
  EdgeKey rev{b, a};
  auto it = edges.find(rev);
  if(it != edges.end()) {
    if(--(it->second) == 0) edges.erase(it);
    return;
  }
  EdgeKey key{a, b};
  ++edges[key];
}

void rebuild_edges(Poly& poly) {
  poly.bb = bbox_of(poly.v);
  poly.H.clear();
  poly.V.clear();
  size_t n = poly.v.size();
  if(n == 0) return;
  for(size_t i = 0; i < n; ++i) {
    const Pt& a = poly.v[i];
    const Pt& b = poly.v[(i + 1) % n];
    if(a.y == b.y) poly.H.push_back({a, b});
    else if(a.x == b.x) poly.V.push_back({a, b});
  }
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

  for(u32 aa_pid = 0; aa_pid < aa_layer.polys.size(); ++aa_pid) {
    const auto& AA = aa_layer.polys[aa_pid];
    std::vector<u32> cands;
    poly_grid.query(AA.bb, cands);
    std::sort(cands.begin(), cands.end());
    cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
    for(u32 poly_pid : cands) {
      if(poly_pid >= poly_layer.polys.size()) continue;
      const auto& P = poly_layer.polys[poly_pid];
      if(!bbox_overlap(AA.bb, P.bb)) continue;
      if(!poly_intersect_manhattan(AA, P)) continue;
      poly_to_aa_cands[poly_pid].push_back(aa_pid);
      aa_to_poly_cands[aa_pid].push_back(poly_pid);
    }
  }
}

void GateCtx::compute_aa_plan(const Layout&, u32 aa_pid) {
  aa_plans[aa_pid].computed = true;
}

void GateCtx::execute_cut(const Layout& L, u32 aa_pid) {
  if(aa_pid >= L.layers[aa_lid].polys.size()) return;
  auto& slices = aa_slices[aa_pid];
  if(slices.ready) return;

  const Poly& AA = L.layers[aa_lid].polys[aa_pid];

  std::vector<i32> xs = {AA.bb.minx, AA.bb.maxx};
  std::vector<i32> ys = {AA.bb.miny, AA.bb.maxy};
  std::vector<Rect> rects;

  for(u32 poly_pid : aa_to_poly_cands[aa_pid]) {
    if(poly_pid >= L.layers[poly_lid].polys.size()) continue;
    const Poly& P = L.layers[poly_lid].polys[poly_pid];
    if(!bbox_overlap(AA.bb, P.bb)) continue;
    if(!poly_intersect_manhattan(AA, P)) continue;
    i32 x0 = std::max(AA.bb.minx, P.bb.minx);
    i32 x1 = std::min(AA.bb.maxx, P.bb.maxx);
    i32 y0 = std::max(AA.bb.miny, P.bb.miny);
    i32 y1 = std::min(AA.bb.maxy, P.bb.maxy);
    if(x0 >= x1 || y0 >= y1) continue;
    Rect r;
    r.x0 = x0; r.x1 = x1; r.y0 = y0; r.y1 = y1;
    r.high = (poly_pid < poly_high.size() && poly_high[poly_pid]);
    rects.push_back(r);
    xs.push_back(x0); xs.push_back(x1);
    ys.push_back(y0); ys.push_back(y1);
  }

  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  std::sort(ys.begin(), ys.end());
  ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

  int nx = std::max<int>(xs.size() - 1, 0);
  int ny = std::max<int>(ys.size() - 1, 0);

  if(nx <= 0 || ny <= 0) {
    Poly piece = AA;
    piece.pid = 0;
    piece.gid = make_slice_gid(aa_pid, 0);
    slices.pieces = {piece};
    slices.adjacency.assign(1, {});
    slices.ready = true;
    return;
  }

  size_t vertical_count = (nx > 1) ? (size_t)(nx - 1) * ny : 0;
  size_t horizontal_count = (ny > 1) ? (size_t)nx * (ny - 1) : 0;
  std::vector<Boundary> vertical(vertical_count);
  std::vector<Boundary> horizontal(horizontal_count);

  auto mark_vertical = [&](int bx, int iy, bool high) {
    if(nx <= 1) return;
    if(bx <= 0 || bx >= (int)xs.size()) return;
    Boundary& b = vertical[(size_t)(bx - 1) * ny + iy];
    b.mark(high);
  };
  auto mark_horizontal = [&](int ix, int by, bool high) {
    if(ny <= 1) return;
    if(by <= 0 || by >= (int)ys.size()) return;
    Boundary& b = horizontal[(size_t)ix * (ny - 1) + (by - 1)];
    b.mark(high);
  };

  for(const auto& r : rects) {
    int ix0 = std::lower_bound(xs.begin(), xs.end(), r.x0) - xs.begin();
    int ix1 = std::lower_bound(xs.begin(), xs.end(), r.x1) - xs.begin();
    int iy0 = std::lower_bound(ys.begin(), ys.end(), r.y0) - ys.begin();
    int iy1 = std::lower_bound(ys.begin(), ys.end(), r.y1) - ys.begin();
    bool high = r.high;

    int bL = std::max(1, ix0);
    int bR = std::min((int)xs.size() - 2, ix1);
    for(int bx = bL; bx <= bR; ++bx) {
      for(int iy = iy0; iy < iy1; ++iy) {
        mark_vertical(bx, iy, high);
      }
    }

    int bB = std::max(1, iy0);
    int bT = std::min((int)ys.size() - 2, iy1);
    for(int by = bB; by <= bT; ++by) {
      for(int ix = ix0; ix < ix1; ++ix) {
        mark_horizontal(ix, by, high);
      }
    }
  }

  auto idx = [&](int ix, int iy) { return ix * ny + iy; };

  std::vector<CellRect> cells;
  cells.reserve((size_t)nx * ny);
  for(int ix = 0; ix < nx; ++ix) {
    for(int iy = 0; iy < ny; ++iy) {
      CellRect cell;
      cell.x0 = xs[ix]; cell.x1 = xs[ix + 1];
      cell.y0 = ys[iy]; cell.y1 = ys[iy + 1];
      cells.push_back(cell);
    }
  }

  std::vector<std::vector<int>> adjacency_cells(cells.size());
  for(int ix = 0; ix < nx; ++ix) {
    for(int iy = 0; iy < ny; ++iy) {
      int u = idx(ix, iy);
      if(nx > 1 && ix + 1 < nx) {
        int bx = ix + 1;
        const Boundary& b = vertical[(size_t)(bx - 1) * ny + iy];
        if(b.open || !b.blocked) {
          int v = idx(ix + 1, iy);
          adjacency_cells[u].push_back(v);
          adjacency_cells[v].push_back(u);
        }
      }
      if(ny > 1 && iy + 1 < ny) {
        int by = iy + 1;
        const Boundary& b = horizontal[(size_t)ix * (ny - 1) + (by - 1)];
        if(b.open || !b.blocked) {
          int v = idx(ix, iy + 1);
          adjacency_cells[u].push_back(v);
          adjacency_cells[v].push_back(u);
        }
      }
    }
  }

  for(auto& nbrs : adjacency_cells) {
    std::sort(nbrs.begin(), nbrs.end());
    nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
  }

  std::vector<int> comp_id(cells.size(), -1);
  std::vector<std::vector<int>> components;
  for(size_t i = 0; i < cells.size(); ++i) {
    if(comp_id[i] != -1) continue;
    std::vector<int> comp;
    std::queue<int> qq;
    comp_id[i] = components.size();
    qq.push((int)i);
    while(!qq.empty()) {
      int u = qq.front(); qq.pop();
      comp.push_back(u);
      for(int v : adjacency_cells[u]) {
        if(comp_id[v] != -1) continue;
        comp_id[v] = comp_id[i];
        qq.push(v);
      }
    }
    components.push_back(std::move(comp));
  }

  std::vector<Poly> comp_polys;
  comp_polys.reserve(components.size());
  for(const auto& comp : components) {
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edges;
    for(int idx_cell : comp) {
      const CellRect& c = cells[idx_cell];
      PointKey a{c.x0, c.y0};
      PointKey b{c.x1, c.y0};
      PointKey c1{c.x1, c.y1};
      PointKey d{c.x0, c.y1};
      add_edge(edges, a, b);
      add_edge(edges, b, c1);
      add_edge(edges, c1, d);
      add_edge(edges, d, a);
    }

    if(edges.empty()) continue;

    std::unordered_map<PointKey, PointKey, PointKeyHash> next_map;
    PointKey start{};
    bool start_set = false;
    for(const auto& kv : edges) {
      if(kv.second <= 0) continue;
      next_map[kv.first.a] = kv.first.b;
      if(!start_set || kv.first.a < start) {
        start = kv.first.a;
        start_set = true;
      }
    }
    if(!start_set) continue;

    std::vector<Pt> poly_pts;
    PointKey cur = start;
    do {
      poly_pts.push_back({cur.x, cur.y});
      auto it = next_map.find(cur);
      if(it == next_map.end()) break;
      cur = it->second;
    } while(cur.x != start.x || cur.y != start.y);

    if(poly_pts.size() < 3) continue;

    Poly poly;
    poly.lid = AA.lid;
    poly.v = std::move(poly_pts);
    rebuild_edges(poly);
    comp_polys.push_back(std::move(poly));
  }

  for(size_t i = 0; i < comp_polys.size(); ++i) {
    comp_polys[i].pid = i;
    comp_polys[i].gid = make_slice_gid(aa_pid, (u32)i);
  }

  slices.pieces = std::move(comp_polys);
  slices.adjacency.assign(slices.pieces.size(), {});
  slices.ready = true;
}

void GateCtx::enqueue_entry_slices_from(u64 u_gid, u32 aa_pid,
                                        std::queue<u64>& q,
                                        std::unordered_set<u64>& vis,
                                        const Layout& L) const {
  auto it = aa_slices.find(aa_pid);
  if(it == aa_slices.end() || !it->second.ready) return;

  const auto& S = it->second;
  const auto& pieces = S.pieces;

  const Poly* source = nullptr;
  if(is_slice_gid(u_gid)) {
    source = slice_from_gid(u_gid);
  } else {
    u32 lid = gid_lid(u_gid);
    u32 pid = gid_pid(u_gid);
    if(lid < L.layers.size() && pid < L.layers[lid].polys.size()) {
      source = &L.layers[lid].polys[pid];
    }
  }
  if(!source) return;

  std::vector<int> entries;
  for(size_t i = 0; i < pieces.size(); ++i) {
    if(!bbox_overlap(pieces[i].bb, source->bb)) continue;
    if(poly_intersect_manhattan(pieces[i], *source)) {
      entries.push_back((int)i);
    }
  }

  if(entries.empty()) {
    long long best = -1;
    int best_idx = -1;
    for(size_t i = 0; i < pieces.size(); ++i) {
      i32 minx = std::max(pieces[i].bb.minx, source->bb.minx);
      i32 maxx = std::min(pieces[i].bb.maxx, source->bb.maxx);
      i32 miny = std::max(pieces[i].bb.miny, source->bb.miny);
      i32 maxy = std::min(pieces[i].bb.maxy, source->bb.maxy);
      long long area = (long long)std::max(0, maxx - minx) * std::max(0, maxy - miny);
      if(area > best) {
        best = area;
        best_idx = (int)i;
      }
    }
    if(best_idx >= 0) entries.push_back(best_idx);
  }

  std::queue<int> iq;
  std::vector<char> seen(pieces.size(), 0);

  for(int idx : entries) {
    u64 gid = make_slice_gid(aa_pid, idx);
    if(vis.insert(gid).second) {
      q.push(gid);
      iq.push(idx);
      seen[idx] = 1;
    }
  }

  while(!iq.empty()) {
    int u = iq.front(); iq.pop();
    if(u >= (int)S.adjacency.size()) continue;
    for(u32 v : S.adjacency[u]) {
      if(v >= pieces.size()) continue;
      if(seen[v]) continue;
      u64 gid = make_slice_gid(aa_pid, v);
      if(vis.insert(gid).second) {
        q.push(gid);
        iq.push((int)v);
        seen[v] = 1;
      }
    }
  }
}

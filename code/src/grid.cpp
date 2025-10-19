#include "grid.hpp"
#include <algorithm>

static void global_bounds(const std::vector<Poly>& polys, i32& minx,i32& miny,i32& maxx,i32& maxy){
  if(polys.empty()){ minx=miny=maxx=maxy=0; return; }
  minx=polys[0].bb.minx; maxx=polys[0].bb.maxx;
  miny=polys[0].bb.miny; maxy=polys[0].bb.maxy;
  for(auto &p: polys){
    minx = std::min(minx, p.bb.minx); maxx = std::max(maxx, p.bb.maxx);
    miny = std::min(miny, p.bb.miny); maxy = std::max(maxy, p.bb.maxy);
  }
}

Grid Grid::build(const std::vector<Poly>& polys, int cell_hint){
  Grid g;
  if(polys.empty()){ return g; }
  // estimate cell size ~ median of bbox widths/heights
  std::vector<int> wh; wh.reserve(polys.size()*2);
  for(auto &p: polys){ wh.push_back(p.bb.maxx - p.bb.minx); wh.push_back(p.bb.maxy - p.bb.miny); }
  std::nth_element(wh.begin(), wh.begin()+wh.size()/2, wh.end());
  int med = std::max(8, wh[wh.size()/2]);
  g.cell = (cell_hint>0? cell_hint : std::max(16, std::min(4096, med)));
  global_bounds(polys, g.minx,g.miny, g.maxx,g.maxy);
  g.nx = ( (g.maxx - g.minx) / g.cell ) + 1;
  g.ny = ( (g.maxy - g.miny) / g.cell ) + 1;
  g.buckets.assign(g.nx*g.ny, {});

  auto cell_id = [&](i32 x, i32 y)->std::pair<int,int>{
    int cx = (x - g.minx) / g.cell;
    int cy = (y - g.miny) / g.cell;
    cx = std::max(0, std::min(g.nx-1, cx));
    cy = std::max(0, std::min(g.ny-1, cy));
    return {cx,cy};
  };

  for(auto &p: polys){
    auto [cx0,cy0] = cell_id(p.bb.minx, p.bb.miny);
    auto [cx1,cy1] = cell_id(p.bb.maxx, p.bb.maxy);
    for(int cy=cy0; cy<=cy1; ++cy){
      for(int cx=cx0; cx<=cx1; ++cx){
        g.buckets[cy*g.nx + cx].push_back(p.pid);
      }
    }
  }
  return g;
}

void Grid::query(const BBox& bb, std::vector<u32>& out) const{
  if(nx==0||ny==0) return;
  auto clamp=[&](int v,int lo,int hi){ return (v<lo?lo:(v>hi?hi:v)); };
  
  // 扩展查询范围以避免边界遗漏
  auto cx0 = clamp((bb.minx - minx)/cell - 1, 0, nx-1);
  auto cy0 = clamp((bb.miny - miny)/cell - 1, 0, ny-1);
  auto cx1 = clamp((bb.maxx - minx)/cell + 1, 0, nx-1);
  auto cy1 = clamp((bb.maxy - miny)/cell + 1, 0, ny-1);
  
  for(int cy=cy0; cy<=cy1; ++cy)
    for(int cx=cx0; cx<=cx1; ++cx){
      const auto& v = buckets[cy*nx + cx];
      out.insert(out.end(), v.begin(), v.end());
    }
}

#pragma once
#include "types.hpp"
#include <vector>

struct Grid {
  int cell=64;          // cell size
  i32 minx=0,miny=0;    // global offset
  i32 maxx=0,maxy=0;    // global extent
  int nx=0, ny=0;       // dims
  // one layer => buckets over pids (dense 0..N-1)
  std::vector<std::vector<u32>> buckets; // size = nx*ny

  static Grid build(const std::vector<Poly>& polys, int cell_hint=-1);

  void query(const BBox& bb, std::vector<u32>& out) const; // MAY produce duplicates
};

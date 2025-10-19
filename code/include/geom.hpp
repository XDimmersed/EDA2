#pragma once
#include "types.hpp"
#include <vector>

BBox bbox_of(const std::vector<Pt>& v);
void  ensure_ccw(std::vector<Pt>& v);
bool  point_in_poly_manhattan(const Pt& p, const Poly& P); // boundary counts as inside
inline bool bbox_overlap(const BBox& a, const BBox& b){
  return !(a.maxx < b.minx || b.maxx < a.minx || a.maxy < b.miny || b.maxy < a.miny);
}
bool segs_intersect_axis_aligned(const Pt& a, const Pt& b, const Pt& c, const Pt& d); // includes endpoints & colinear overlap
bool poly_intersect_manhattan(const Poly& A, const Poly& B); // area overlap OR boundary (point/segment) intersection

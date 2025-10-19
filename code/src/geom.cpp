#include "geom.hpp"
#include <algorithm>

static long long signed_area2(const std::vector<Pt>& v){
  long long s=0; int n=(int)v.size();
  for(int i=0;i<n;i++){ int j=(i+1)%n; s += 1LL*v[i].x*v[j].y - 1LL*v[j].x*v[i].y; }
  return s;
}

BBox bbox_of(const std::vector<Pt>& v){
  BBox b; if(v.empty()) return b;
  b.minx=b.maxx=v[0].x; b.miny=b.maxy=v[0].y;
  for(auto &p: v){
    b.minx = std::min(b.minx, p.x);
    b.maxx = std::max(b.maxx, p.x);
    b.miny = std::min(b.miny, p.y);
    b.maxy = std::max(b.maxy, p.y);
  }
  return b;
}

void ensure_ccw(std::vector<Pt>& v){
  if(v.size()<3) return;
  if(signed_area2(v) < 0) std::reverse(v.begin(), v.end());
}

static inline bool on_segment_colinear(i32 a,i32 b,i32 c){
  return std::min(a,b)<=c && c<=std::max(a,b);
}

bool segs_intersect_axis_aligned(const Pt& a, const Pt& b, const Pt& c, const Pt& d){
  // 允许端点接触与共线重合；适配 H×V、V×H 以及同向共线
  bool abH = (a.y==b.y), cdH = (c.y==d.y);
  if(abH && !cdH){
    // A水平，C垂直
    i32 y=a.y; i32 x=c.x;
    if(on_segment_colinear(a.x,b.x,x) && on_segment_colinear(c.y,d.y,y)) return true;
  }else if(!abH && cdH){
    // A垂直，C水平
    i32 x=a.x; i32 y=c.y;
    if(on_segment_colinear(c.x,d.x,x) && on_segment_colinear(a.y,b.y,y)) return true;
  }else if(abH && cdH){
    // 同为水平且共线时的重叠
    if(a.y==c.y){
      if(std::max(std::min(a.x,b.x), std::min(c.x,d.x)) <=
         std::min(std::max(a.x,b.x), std::max(c.x,d.x))) return true;
    }
  }else{
    // 同为垂直且共线时的重叠
    if(a.x==c.x){
      if(std::max(std::min(a.y,b.y), std::min(c.y,d.y)) <=
         std::min(std::max(a.y,b.y), std::max(c.y,d.y))) return true;
    }
  }
  return false;
}

static bool point_on_edge(const Pt& p, const Pt& a, const Pt& b){
  if(a.x==b.x){ // 垂直边
    if(p.x!=a.x) return false;
    return on_segment_colinear(a.y,b.y,p.y);
  }else if(a.y==b.y){ // 水平边
    if(p.y!=a.y) return false;
    return on_segment_colinear(a.x,b.x,p.x);
  }
  return false; // 理论上不会出现（曼哈顿）
}

bool point_in_poly_manhattan(const Pt& p, const Poly& P){
  if(p.x < P.bb.minx || p.x > P.bb.maxx || p.y < P.bb.miny || p.y > P.bb.maxy) return false;
  const auto& v = P.v; int n=(int)v.size();
  for(int i=0;i<n;i++){
    int j=(i+1)%n;
    if(point_on_edge(p, v[i], v[j])) return true;
  }
  bool inside=false;
  for(int i=0;i<n;i++){
    const Pt& a = v[i];
    const Pt& b = v[(i+1)%n];
    if(a.x!=b.x) continue; // 仅处理垂直边
    i32 miny = std::min(a.y,b.y);
    i32 maxy = std::max(a.y,b.y);
    if(p.y < miny || p.y >= maxy) continue; // 下闭上开，避免顶点重复计数
    if(a.x <= p.x) continue;                // 射线 (+x) 右侧交点
    inside = !inside;
  }
  return inside;
}

bool poly_intersect_manhattan(const Poly& A, const Poly& B){
  if(!bbox_overlap(A.bb, B.bb)) return false;
  
  // 1. 边相交检查：H×V / V×H
  for(auto &eA : A.H){
    for(auto &eB : B.V){
      if(segs_intersect_axis_aligned(eA.first, eA.second, eB.first, eB.second)) return true;
    }
  }
  for(auto &eA : A.V){
    for(auto &eB : B.H){
      if(segs_intersect_axis_aligned(eA.first, eA.second, eB.first, eB.second)) return true;
    }
  }
  
  // 2. 顶点包含检查（包括边界）
  for(const auto& pt : A.v) {
    if(point_in_poly_manhattan(pt, B)) return true;
  }
  for(const auto& pt : B.v) {
    if(point_in_poly_manhattan(pt, A)) return true;
  }
  
  // 3. 边界接触检查：A的顶点在B的边上，或B的顶点在A的边上
  for(const auto& pt : A.v) {
    for(size_t i = 0; i < B.v.size(); ++i) {
      const Pt& b1 = B.v[i];
      const Pt& b2 = B.v[(i+1) % B.v.size()];
      if(point_on_edge(pt, b1, b2)) return true;
    }
  }
  for(const auto& pt : B.v) {
    for(size_t i = 0; i < A.v.size(); ++i) {
      const Pt& a1 = A.v[i];
      const Pt& a2 = A.v[(i+1) % A.v.size()];
      if(point_on_edge(pt, a1, a2)) return true;
    }
  }
  
  // 4. 完全覆盖兜底：A 内一点在 B 内或反之
  if(!A.v.empty() && point_in_poly_manhattan(A.v[0], B)) return true;
  if(!B.v.empty() && point_in_poly_manhattan(B.v[0], A)) return true;
  
  return false;
}

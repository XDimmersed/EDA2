#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

using i32 = int32_t; using u32 = uint32_t; using u64 = uint64_t;

struct Pt { i32 x=0, y=0; };
struct BBox { i32 minx=0, miny=0, maxx=0, maxy=0; };

struct Poly {
  u32 lid=0;                 // dense layer id
  u32 pid=0;                 // dense poly id within layer
  u64 gid=0;                 // ((u64)lid<<32) | pid
  BBox bb;
  std::vector<Pt> v;         // CCW vertices
  std::vector<std::pair<Pt,Pt>> H, V; // horizontal & vertical edges
};

struct Layer {
  u32 lid=0;
  std::string name;
  std::vector<Poly> polys;
};

struct Rule {
  struct Seed { std::string layer; Pt p; };
  std::vector<Seed> seeds;                       // 1~2 seeds
  std::vector<std::vector<std::string>> vias;    // each is a chain, adjacent-only connectivity
  bool has_gate=false;
  std::string poly_layer_name;                   // Gate: first = poly
  std::string aa_layer_name;                     // Gate: second = AA
};

struct Layout {
  std::vector<Layer> layers;
  std::unordered_map<std::string,u32> name2lid;
};

struct Config {
  int threads = 1;
  int grid_cell_hint = -1;   // auto if <0
  bool use_hot_rtree = false; // placeholder hook
};

// helpers
inline u64 make_gid(u32 lid, u32 pid){ return ( (u64)lid<<32 ) | (u64)pid; }
inline u32 gid_lid(u64 gid){ return (u32)(gid>>32); }
inline u32 gid_pid(u64 gid){ return (u32)(gid & 0xffffffffULL); }

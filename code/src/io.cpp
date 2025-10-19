#include "io.hpp"
#include <fstream>
#include <unordered_map>
#include <algorithm>

static const Poly* resolve_poly(const Layout& L, const GateCtx* gate, u64 gid){
  u32 lid = gid_lid(gid);
  u32 pid = gid_pid(gid);
  if(gate && gate->is_slice_gid(gid)){
    return gate->slice_from_gid(gid);
  }
  if(lid>=L.layers.size()) return nullptr;
  const auto& layer = L.layers[lid];
  if(pid>=layer.polys.size()) return nullptr;
  return &layer.polys[pid];
}

void write_result(const Layout& L, const std::vector<u64>& reachable, const GateCtx* gate, const std::string& out_path){
  std::unordered_map<u32, std::vector<u64>> layer_gids;
  layer_gids.reserve(L.layers.size());
  for(u64 gid : reachable){
    const Poly* poly = resolve_poly(L, gate, gid);
    if(!poly || poly->v.size()<3) continue;
    layer_gids[poly->lid].push_back(gid);
  }
  for(auto &kv : layer_gids){
    auto &vec = kv.second;
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  }

  std::ofstream fout(out_path);
  for(const auto& layer : L.layers){
    auto it = layer_gids.find(layer.lid);
    if(it==layer_gids.end() || it->second.empty()) continue;
    fout << layer.name << "\n";
    for(u64 gid : it->second){
      const Poly* poly = resolve_poly(L, gate, gid);
      if(!poly || poly->v.size()<3) continue;
      for(size_t i=0;i<poly->v.size();++i){
        fout << "(" << poly->v[i].x << "," << poly->v[i].y << ")";
        if(i+1<poly->v.size()) fout << ",";
      }
      fout << "\n";
    }
  }
}

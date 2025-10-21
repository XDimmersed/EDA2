#include "parser.hpp"
#include "fdec.hpp"
#include "io.hpp"
#include "equiv_class.hpp"
#include <iostream>

int main(int argc, char** argv){
  std::string layout_path, rule_path, out_path;
  int threads=1;

  for(int i=1;i<argc;i++){
    std::string a=argv[i];
    auto need=[&](const char* k){ return a==k && i+1<argc; };
    if(need(std::string("-layout").c_str())) layout_path=argv[++i];
    else if(need(std::string("-rule").c_str())) rule_path=argv[++i];
    else if(need(std::string("-output").c_str())) out_path=argv[++i];
    else if(need(std::string("-thread").c_str())) threads=std::stoi(argv[++i]);
  }
  if(layout_path.empty()||rule_path.empty()||out_path.empty()){
    std::cerr<<"Usage: trace -layout <layout.txt> -rule <rule.txt> [-thread n] -output <res.txt>\n";
    return 1;
  }

  Rule   R = parse_rule_file(rule_path);
  Layout L = parse_layout_file(layout_path);
  Config cfg; cfg.threads = threads;

  // std::cerr << "Parsed " << L.layers.size() << " layers, " << R.seeds.size() << " seeds\n";
  // std::cerr << "Has gate: " << (R.has_gate ? "yes" : "no") << "\n";

  // 诊断输入数据的几何重复（在任何BFS之前）
  analyze_input_duplicates(L);

  // 构建等价类（预扫描）
  EquivClassManager eqv_mgr;
  u32 aa_lid = UINT32_MAX, poly_lid = UINT32_MAX;
  if(R.has_gate) {
    auto it1 = L.name2lid.find(R.poly_layer_name);
    auto it2 = L.name2lid.find(R.aa_layer_name);
    if(it1 != L.name2lid.end()) poly_lid = it1->second;
    if(it2 != L.name2lid.end()) aa_lid = it2->second;
  }
  eqv_mgr.build_classes(&L, aa_lid, poly_lid);

  FDEC eng; eng.build(L, R, cfg);

  std::vector<u64> reachable;
  if(!R.has_gate || R.seeds.size()<=1){
    // 第一/第二问：取所有 seed 作为多源
    std::vector<std::pair<u32,u32>> seeds;
    for(auto &s : R.seeds){
      auto it = L.name2lid.find(s.layer);
      if(it==L.name2lid.end()) continue;
      u32 lid = it->second;
      u32 pid = eng.locate_seed_pid(lid, s.p);
      if(pid!=UINT32_MAX) seeds.push_back({lid,pid});
    }
    reachable = eng.trace_no_gate(seeds);
  }else{
    // 第三问
    reachable = eng.trace_with_gate();
    std::cerr << "Found " << reachable.size() << " reachable polygons (with gate)\n";
  }

  write_result(L, reachable, eng.gate_result(), out_path);
  return 0;
}

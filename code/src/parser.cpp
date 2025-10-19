#include "parser.hpp"
#include "geom.hpp"
#include <fstream>
#include <sstream>
#include <cctype>
#include <stdexcept>
#include <unordered_map>
#include <iostream>

static std::string trim(const std::string& s){
  size_t i=0,j=s.size();
  while(i<j && std::isspace((unsigned char)s[i])) ++i;
  while(j>i && std::isspace((unsigned char)s[j-1])) --j;
  return s.substr(i,j-i);
}

static bool parse_point(const std::string& tok, Pt& p){
  // format: (x,y)
  if(tok.size()<5) return false;
  if(tok.front()!='(' || tok.back()!=')') return false;
  auto inner = tok.substr(1, tok.size()-2);
  auto comma = inner.find(',');
  if(comma==std::string::npos) return false;
  p.x = std::stoi(inner.substr(0, comma));
  p.y = std::stoi(inner.substr(comma+1));
  return true;
}

Rule parse_rule_file(const std::string& rule_path){
  std::ifstream fin(rule_path);
  if(!fin) throw std::runtime_error("Cannot open rule file: "+rule_path);
  Rule R;
  std::string line;
  enum Sect{NONE, STARTPOS, VIA, GATE} sect=NONE;
  while(std::getline(fin,line)){
    line = trim(line);
    if(line.empty()) continue;
    if(line=="StartPos"){ sect=STARTPOS; continue; }
    if(line=="Via"){ sect=VIA; continue; }
    if(line=="Gate"){ sect=GATE; continue; }
    if(sect==STARTPOS){
      std::istringstream ss(line);
      std::string lname; 
      ss>>lname;
      // Parse coordinate in format (x,y)
      std::string coord_str;
      ss>>coord_str;
      Pt p;
      if(parse_point(coord_str, p)){
        R.seeds.push_back({lname, p});
      }
    }else if(sect==VIA){
      std::istringstream ss(line);
      std::string name; std::vector<std::string> chain;
      while(ss>>name) chain.push_back(name);
      if(!chain.empty()) R.vias.push_back(chain);
    }else if(sect==GATE){
      std::istringstream ss(line);
      ss>>R.poly_layer_name>>R.aa_layer_name;
      R.has_gate = (!R.poly_layer_name.empty() && !R.aa_layer_name.empty());
    }
  }
  return R;
}

Layout parse_layout_file(const std::string& layout_path){
  std::ifstream fin(layout_path);
  if(!fin) throw std::runtime_error("Cannot open layout file: "+layout_path);

  Layout L;
  std::string line, cur_layer;
  Layer* layer=nullptr;
  std::unordered_map<std::string,u32>& name2lid = L.name2lid;
  int line_num = 0;
  
  while(std::getline(fin,line)){
    line_num++;
    line = trim(line);
    if(line.empty()) continue;
    
    // Layer name line: token without parentheses
    if(line.find('(')==std::string::npos){
      cur_layer = line;
      if(!name2lid.count(cur_layer)){
        u32 lid = (u32)L.layers.size();
        name2lid[cur_layer]=lid;
        L.layers.push_back(Layer{lid,cur_layer,{}});
      }
      layer = &L.layers[name2lid[cur_layer]];
      continue;
    }
    
    // polygon line: "(x,y),(x,y),..."
    if(!layer) continue;
    std::vector<Pt> v;
    size_t pos = 0;
    while(pos < line.length()){
      size_t start = line.find('(', pos);
      if(start == std::string::npos) break;
      size_t end = line.find(')', start);
      if(end == std::string::npos) break;
      std::string tok = line.substr(start, end - start + 1);
      Pt p; 
      if(parse_point(tok, p)) v.push_back(p);
      pos = end + 1;
    }
    if(v.size()<3) continue;
    ensure_ccw(v);
    Poly P;
    P.lid = layer->lid;
    P.pid = (u32)layer->polys.size();
    P.gid = make_gid(P.lid, P.pid);
    P.v = std::move(v);
    P.bb = bbox_of(P.v);
    // build H/V edge lists
    int n=(int)P.v.size();
    for(int i=0;i<n;i++){
      Pt a=P.v[i], b=P.v[(i+1)%n];
      if(a.y==b.y) P.H.push_back({a,b});
      else if(a.x==b.x) P.V.push_back({a,b});
      else {
        // 按题意应为曼哈顿；但若存在非轴对齐边，为健壮性仍收录到 H/V 最近方向（可在后续替换为严格检查）
        if(std::abs(a.y-b.y) < std::abs(a.x-b.x)) P.H.push_back({a,b}); else P.V.push_back({a,b});
      }
    }
    layer->polys.push_back(std::move(P));
  }
  return L;
}

#pragma once
#include "types.hpp"
#include <map>
#include <set>
#include <iostream>
#include <iomanip>

// 度量仪表盘：跟踪整个算法流程的关键指标

struct Metrics {
  // Phase 1统计
  std::map<u32, int> phase1_visited_by_layer;
  std::map<u32, int> phase1_high_by_layer;
  int high_poly_intersect_aa = 0;  // 高电平Poly与AA真实相交数
  int high_poly_penetrate_aa = 0;  // 高电平Poly贯穿AA数
  
  // Phase 2统计
  std::set<u32> encountered_aa;  // 遇到的AA
  std::set<u32> cut_aa;          // 被切的AA
  int total_aa = 0;
  
  // 切割线统计（合并前）
  int raw_vlines_high = 0;
  int raw_vlines_low = 0;
  int raw_hlines_high = 0;
  int raw_hlines_low = 0;
  
  // 切割线统计（合并后）
  int merged_vlines_high = 0;
  int merged_vlines_low = 0;
  int merged_hlines_high = 0;
  int merged_hlines_low = 0;
  
  // 片段统计
  std::map<u32, int> aa_pieces_count;       // AA pid -> 片段数
  std::map<u32, bool> aa_has_high_adj;      // AA pid -> 是否有高电平邻接
  std::vector<std::pair<u32, int>> isolated_slices_distances;  // 孤立片段到最近高电平线的距离
  
  // 输出统计
  std::map<std::string, int> output_slices;    // 层名 -> 切片数量
  std::map<std::string, int> output_originals; // 层名 -> 原始多边形数量
  std::map<std::string, int> output_dropped;   // 层名 -> 去重数量
  
  void print_phase1_summary() const {
    std::cerr << "\n========================================\n";
    std::cerr << "  Phase 1 度量仪表盘\n";
    std::cerr << "========================================\n";
    
    std::cerr << "\n[Phase 1 按层统计]\n";
    int total_visited = 0, total_high = 0;
    for(const auto& [lid, count] : phase1_visited_by_layer) {
      total_visited += count;
      int high = phase1_high_by_layer.count(lid) ? phase1_high_by_layer.at(lid) : 0;
      total_high += high;
      if(count > 0) {
        std::cerr << "  Layer " << lid << ": visited=" << count 
                  << ", poly_high=" << high 
                  << " (" << (count > 0 ? high*100/count : 0) << "%)\n";
      }
    }
    std::cerr << "  Total: visited=" << total_visited 
              << ", poly_high=" << total_high << "\n";
    
    std::cerr << "\n[高电平Poly与AA关系]\n";
    std::cerr << "  高电平Poly与AA真实相交: " << high_poly_intersect_aa << "\n";
    std::cerr << "  高电平Poly贯穿AA: " << high_poly_penetrate_aa << "\n";
    if(high_poly_intersect_aa > 0) {
      std::cerr << "  贯穿率: " << (high_poly_penetrate_aa*100/high_poly_intersect_aa) << "%\n";
    }
  }
  
  void print_phase2_summary() const {
    std::cerr << "\n========================================\n";
    std::cerr << "  Phase 2 度量仪表盘\n";
    std::cerr << "========================================\n";
    
    std::cerr << "\n[AA覆盖统计]\n";
    std::cerr << "  全部AA数: " << total_aa << "\n";
    std::cerr << "  遇到的AA数: " << encountered_aa.size() 
              << " (" << (total_aa > 0 ? encountered_aa.size()*100/total_aa : 0) << "%)\n";
    std::cerr << "  被切的AA数: " << cut_aa.size() 
              << " (" << (encountered_aa.size() > 0 ? cut_aa.size()*100/encountered_aa.size() : 0) << "% of encountered)\n";
    std::cerr << "  未遇到的AA: " << (total_aa - encountered_aa.size()) << "\n";
    std::cerr << "  遇到但未切: " << (encountered_aa.size() - cut_aa.size()) << "\n";
    
    std::cerr << "\n[切割线统计 - 合并前]\n";
    int raw_v_total = raw_vlines_high + raw_vlines_low;
    int raw_h_total = raw_hlines_high + raw_hlines_low;
    std::cerr << "  垂直线: " << raw_v_total 
              << " (高=" << raw_vlines_high << ", 低=" << raw_vlines_low;
    if(raw_v_total > 0) {
      std::cerr << ", 高占比=" << (raw_vlines_high*100/raw_v_total) << "%";
    }
    std::cerr << ")\n";
    std::cerr << "  水平线: " << raw_h_total 
              << " (高=" << raw_hlines_high << ", 低=" << raw_hlines_low;
    if(raw_h_total > 0) {
      std::cerr << ", 高占比=" << (raw_hlines_high*100/raw_h_total) << "%";
    }
    std::cerr << ")\n";
    
    std::cerr << "\n[切割线统计 - 合并后]\n";
    int merged_v_total = merged_vlines_high + merged_vlines_low;
    int merged_h_total = merged_hlines_high + merged_hlines_low;
    std::cerr << "  垂直线: " << merged_v_total 
              << " (高=" << merged_vlines_high << ", 低=" << merged_vlines_low;
    if(merged_v_total > 0) {
      std::cerr << ", 高占比=" << (merged_vlines_high*100/merged_v_total) << "%";
    }
    std::cerr << ")\n";
    std::cerr << "  水平线: " << merged_h_total 
              << " (高=" << merged_hlines_high << ", 低=" << merged_hlines_low;
    if(merged_h_total > 0) {
      std::cerr << ", 高占比=" << (merged_hlines_high*100/merged_h_total) << "%";
    }
    std::cerr << ")\n";
    
    if(raw_v_total > 0 || raw_h_total > 0) {
      std::cerr << "  合并率: " << (100 - (merged_v_total + merged_h_total)*100/(raw_v_total + raw_h_total)) << "%\n";
    }
  }
  
  void print_slices_summary() const {
    std::cerr << "\n[片段统计]\n";
    
    // 统计片段分布
    std::map<int, int> pieces_distribution;  // 片段数 -> AA数量
    int total_has_high_adj = 0;
    for(const auto& [aa_pid, pieces] : aa_pieces_count) {
      pieces_distribution[pieces]++;
      if(aa_has_high_adj.count(aa_pid) && aa_has_high_adj.at(aa_pid)) {
        total_has_high_adj++;
      }
    }
    
    std::cerr << "  AA片段分布:\n";
    for(const auto& [pieces, count] : pieces_distribution) {
      std::cerr << "    " << pieces << "片: " << count << " AA\n";
    }
    
    std::cerr << "  有高电平邻接的AA: " << total_has_high_adj 
              << " / " << aa_pieces_count.size();
    if(aa_pieces_count.size() > 0) {
      std::cerr << " (" << (total_has_high_adj*100/aa_pieces_count.size()) << "%)";
    }
    std::cerr << "\n";
    
    if(!isolated_slices_distances.empty()) {
      std::cerr << "\n  孤立片段分析 (前10个):\n";
      int show_count = std::min(10, (int)isolated_slices_distances.size());
      for(int i = 0; i < show_count; ++i) {
        auto [slice_gid, dist] = isolated_slices_distances[i];
        std::cerr << "    片段" << slice_gid << ": 最近高电平线距离=" << dist << "\n";
      }
    }
  }
  
  void print_output_summary() const {
    std::cerr << "\n========================================\n";
    std::cerr << "  输出度量仪表盘\n";
    std::cerr << "========================================\n";
    
    std::cerr << "\n[各层输出统计]\n";
    std::cerr << "  层名     切片    原始    去重    总输出\n";
    std::cerr << "  ----------------------------------------\n";
    
    std::set<std::string> all_layers;
    for(const auto& [layer, _] : output_slices) all_layers.insert(layer);
    for(const auto& [layer, _] : output_originals) all_layers.insert(layer);
    
    int total_slices = 0, total_originals = 0, total_dropped = 0;
    for(const auto& layer : all_layers) {
      int slices = output_slices.count(layer) ? output_slices.at(layer) : 0;
      int originals = output_originals.count(layer) ? output_originals.at(layer) : 0;
      int dropped = output_dropped.count(layer) ? output_dropped.at(layer) : 0;
      int total = slices + originals;
      
      std::cerr << "  " << std::left << std::setw(8) << layer 
                << std::right << std::setw(6) << slices
                << std::setw(8) << originals
                << std::setw(8) << dropped
                << std::setw(10) << total << "\n";
      
      total_slices += slices;
      total_originals += originals;
      total_dropped += dropped;
    }
    
    std::cerr << "  ----------------------------------------\n";
    std::cerr << "  Total   " 
              << std::setw(6) << total_slices
              << std::setw(8) << total_originals
              << std::setw(8) << total_dropped
              << std::setw(10) << (total_slices + total_originals) << "\n";
  }
};

// 全局度量对象
extern Metrics g_metrics;


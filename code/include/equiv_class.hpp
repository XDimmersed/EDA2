#pragma once
#include "types.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>

// 去重策略
enum DedupMode {
  NONE = 0,    // 不去重，输出所有成员
  LIGHT = 1,   // 轻度去重（仅GID去重）
  STRICT = 2   // 严格去重（几何去重）
};

// 几何键（用于哈希）
struct PolyKey {
  u64 h1, h2;
  bool operator==(const PolyKey& o) const { return h1 == o.h1 && h2 == o.h2; }
  bool operator<(const PolyKey& o) const {
    if(h1 != o.h1) return h1 < o.h1;
    return h2 < o.h2;
  }
};

struct PolyHasher {
  size_t operator()(const PolyKey& k) const {
    return std::hash<u64>()(k.h1) ^ (std::hash<u64>()(k.h2) << 1);
  }
};

// 等价类信息
struct ClassInfo {
  u32 lid;                     // 层ID
  std::vector<u32> members;    // 等价类成员pid列表
  u32 rep;                     // 代表pid（通常是第一个成员）
  PolyKey key;                 // 几何键
  
  ClassInfo() : lid(0), rep(0) {}
  ClassInfo(u32 l, std::vector<u32> m, u32 r, PolyKey k) 
    : lid(l), members(std::move(m)), rep(r), key(k) {}
};

// 等价类管理器
class EquivClassManager {
public:
  EquivClassManager() : L(nullptr) {}
  
  // 预扫描：建立等价类
  void build_classes(const Layout* layout, u32 aa_layer_id, u32 poly_layer_id);
  
  // 查询：给定gid，返回对应的类ID（如果该层启用了等价类）
  // 返回值：如果是等价类节点，返回true并填充cid；否则返回false
  bool try_get_class_id(u64 gid, u32& out_cid) const;
  
  // 查询：给定层和cid，返回类信息
  const ClassInfo& get_class(u32 lid, u32 cid) const;
  
  // 查询：给定层和cid，返回所有成员pid
  const std::vector<u32>& get_members(u32 lid, u32 cid) const;
  
  // 查询：给定层和cid，返回代表pid
  u32 get_rep_pid(u32 lid, u32 cid) const;
  
  // 查询：该层是否启用了等价类
  bool is_class_layer(u32 lid) const;
  
  // 查询：该层的去重模式
  DedupMode get_layer_mode(u32 lid) const;
  
  // 统计和日志
  void print_statistics() const;
  
  // 获取层的类数量
  size_t get_class_count(u32 lid) const;
  
private:
  const Layout* L;
  
  // layer -> [ClassInfo]
  std::vector<std::vector<ClassInfo>> layer_classes;
  
  // gid -> cid（只对输入层且启用等价类的填充）
  std::unordered_map<u64, u32> gid_to_cid;
  
  // layer -> mode
  std::vector<DedupMode> layer_modes;
  
  // 标记哪些层启用了等价类
  std::unordered_set<u32> class_enabled_layers;
  
  // AA层和Poly层的lid
  u32 aa_lid = UINT32_MAX;
  u32 poly_lid = UINT32_MAX;
  
  // 辅助函数：计算多边形的规范化几何键
  PolyKey compute_canonical_key(const std::vector<Pt>& v) const;
};


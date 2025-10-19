# gate_fixed.cpp Review

## Summary of observed behaviour
- 用户回归结果显示：AA 切割数量从 6 降为 0，L1 层多边形直接清零，同时其它层也大幅缩减，说明新的切割判据几乎禁止了所有 AA-Poly 贯穿的识别。

## Root-cause findings
1. **贯穿判据过度收紧**  
   `is_full_penetration` 要求交集同时覆盖 AA 的两条对向边，还要检测 poly 的 bbox 必须越过这两条边之外。一旦 poly 在某一侧刚好齐平或者只有局部触达，整条贯穿就会被判为“非贯穿”，导致所有切割线都被过滤掉。【F:code/src/gate.cpp†L24-L51】

2. **半格裁剪仍未生成交点**  
   `clip_axis_aligned_halfgrid` 在遇到奇数坐标时依然选择“交点不输出”，寄希望于端点来闭合轮廓；这会把穿越半格切线的垂直/水平边直接截断，生成的切片缺边、不闭合，后续 bbox 映射自然失败。【F:code/src/gate.cpp†L58-L123】

3. **切片未做规范化**  
   切割后只写入 `v/bb/lid` 就直接进入导通阶段，没有执行 `ensure_ccw`、`build_edges` 等解析阶段的标准化流程，从而导致 `poly_intersect_manhattan` 针对入口片的精判继续不可靠。【F:code/src/gate.cpp†L397-L462】【F:code/src/gate.cpp†L640-L688】

4. **仍然依赖网格单元一一映射**  
   即使前面尝试处理半格切线，导通边的重建依旧要求片段 bbox 完全等于一个 `[X[i],X[i+1]]×[Y[j],Y[j+1]]` 单元，稍微产生 L 形或跨多个格子的切片就会落入 `unmapped_count`，高电平导通关系无法建立。【F:code/src/gate.cpp†L464-L582】

## Recommendations
- 重新实现贯穿检测：在 AA 内部沿法向取样（如中心扫描），通过区间覆盖判断 poly 是否在 AA 的法向上形成贯穿，而不是仅看 bbox 是否出界。
- 把 `clip_axis_aligned_halfgrid` 替换为标准 Sutherland–Hodgman 变体：无论整数还是半格切线都应显式输出交点，再在末尾去重首尾重复点以确保多边形闭合。
- 对每个新切片调用与解析阶段一致的规范化例程（确保逆时针、重建水平/垂直边、刷新 `pid/gid`），让几何运算与入口判定有正确的输入。
- 导通边应在切割过程中直接继承父片段并把“高电平切线的两个儿子”互联，彻底摆脱对完美网格的依赖；必要时在切片集合上再做一次几何邻接兜底。

按照以上方向回调，可以恢复 AA 层的切割数量，并让第三问的 BFS 拓扑与高电平门结构匹配。

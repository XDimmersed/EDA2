# 第三问模块回归测试

## 测试目的
对 Gate 模块（第三问）进行快速回归测试，确认 `poly_high` 标记、AA 切割以及 Lazy-Cut 连通逻辑是否按题面工作。

## 测试环境
- 可执行文件：`build/trace`
- 构建方式：`cmake .. && make -j4`
- 运行平台：Linux (容器环境)

## 测试 1：垂直 Poly 贯穿 AA

### 输入数据
`layout.txt`
```text
POLY
(4,-2),(6,-2),(6,6),(4,6)
AA
(0,0),(10,0),(10,4),(0,4)
METAL1
(0,0),(3,0),(3,4),(0,4)
METAL1
(7,0),(10,0),(10,4),(7,4)
```

`rule_high.txt`
```text
StartPos
POLY (5,0)
METAL1 (1,2)
Via
POLY AA
AA METAL1
Gate
POLY AA
```

`rule_low.txt`
```text
StartPos
METAL1 (1,2)
METAL1 (1,2)
Via
AA METAL1
Gate
POLY AA
```

### 执行命令
```bash
./build/trace -layout layout.txt -rule rule_high.txt -output out_high.txt
./build/trace -layout layout.txt -rule rule_low.txt -output out_low.txt
```

### 观察
- 高电平（`rule_high.txt`）时，Phase 1 日志显示正确标记 1 个 Poly，AA 被切成两块，但切割方向是**水平**（AA 输出的两块上下排列）。【4ecf15†L1-L8】【cb59b5†L1-L8】
- 低电平（`rule_low.txt`）时，尽管 Phase 1 日志显示 `poly_high` 为 0，AA 仍被切成同样的两块，上下两块都与左右 Metal1 接触，导致右侧 Metal1 仍然可达，Gate 逻辑失效。【62302b†L1-L8】【b3036e†L1-L7】

## 测试 2：水平 Poly 贯穿 AA

### 输入数据
`layout_horizontal.txt`
```text
POLY
(0,1),(10,1),(10,3),(0,3)
AA
(0,0),(10,0),(10,4),(0,4)
METAL1
(0,0),(10,0),(10,1),(0,1)
METAL1
(0,3),(10,3),(10,4),(0,4)
```

`rule_horizontal_high.txt`
```text
StartPos
POLY (5,2)
METAL1 (5,0)
Via
POLY AA
AA METAL1
Gate
POLY AA
```

`rule_horizontal_low.txt`
```text
StartPos
METAL1 (5,0)
METAL1 (5,0)
Via
AA METAL1
Gate
POLY AA
```

### 执行命令
```bash
./build/trace -layout layout_horizontal.txt -rule rule_horizontal_high.txt -output out_horizontal_high.txt
./build/trace -layout layout_horizontal.txt -rule rule_horizontal_low.txt -output out_horizontal_low.txt
```

### 观察
- 高电平时，AA 被切成左右两块（垂直切割）。【1f1a68†L1-L8】【dd5945†L1-L8】
- 低电平时，AA 仍被切成相同的左右两块，导致上下 Metal1 依旧互通，Gate 未能断开源漏。【a02904†L1-L8】【a9b96a†L1-L7】

## 初次测试结论（修正前）
`GateCtx::compute_aa_plan` 在判断 Poly 贯穿方向时，把“触底+触顶”映射为**水平切割**，把“触左+触右”映射为**垂直切割**，与物理语义相反，导致无论 Poly 是否高电平，AA 的切片都与错误的两侧接触，从而无法隔离源漏。【F:src/gate_new_gpt.cpp†L274-L309】

建议修正贯穿方向判断，使垂直贯穿产生垂直切线、水平贯穿产生水平切线，并补充针对低电平 Poly 的断开用例。

---

## 修正后的复测（2024-10-20）

### 复测目的
在“Fix AA cut orientation for gate planning” 提交之后，再次复现上面的两组用例，确认 AA 层是否会被切割、导通状态是否符合预期。

### 执行命令
```bash
./build/trace -layout /tmp/layout.txt -rule /tmp/rule_high.txt -output /tmp/out_high.txt
./build/trace -layout /tmp/layout.txt -rule /tmp/rule_low.txt -output /tmp/out_low.txt
./build/trace -layout /tmp/layout_horizontal.txt -rule /tmp/rule_horizontal_high.txt -output /tmp/out_horizontal_high.txt
./build/trace -layout /tmp/layout_horizontal.txt -rule /tmp/rule_horizontal_low.txt -output /tmp/out_horizontal_low.txt
```

### 复测观察
- 垂直贯穿用例（`layout.txt`）：Phase 2 现在正确统计到“cut 1 AA polygons”，生成两块纵向切片；高电平规则访问到两块 AA，低电平规则仅访问到左侧入口片，输出文件也分别呈现整块/单侧 AA。【bb9adb†L1-L7】【3deed4†L1-L5】【bceda0†L1-L9】【fd94cb†L1-L5】
- 水平贯穿用例（`layout_horizontal.txt`）：日志同样显示切割 1 块 AA，AA 输出被分成上下两片；低电平规则只保留下方片段与起点导通，上方被阻断。【2904cb†L1-L8】【f06522†L1-L8】【31e1bd†L1-L9】【8f9f06†L1-L5】

### 复测结论
`compute_aa_plan` 在上述场景中能生成正确的切割线，AA 层被分成两块；高电平 Poly 通过导通边联通双片，低电平 Poly 阻断另一侧。此前“cut 0”日志属于统计错误，现已在 Phase 2 末尾改为按切片缓存汇总，可准确反映 AA 切割数量。

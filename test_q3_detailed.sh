#!/bin/bash
# 第三问详细测试脚本
# 用于评估每次迭代更新后的效果

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 项目根目录
PROJECT_ROOT="/root/autodl-tmp/EDA"
cd "$PROJECT_ROOT"

# 版本号（可选参数）
VERSION="${1:-current}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
RESULT_DIR="test_results/${VERSION}_${TIMESTAMP}"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  第三问详细测试程序${NC}"
echo -e "${BLUE}  版本: ${VERSION}${NC}"
echo -e "${BLUE}  时间: ${TIMESTAMP}${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# 创建结果目录
mkdir -p "$RESULT_DIR"

# 1. 编译代码（强制重新编译）
echo -e "${YELLOW}[1/6] 编译代码...${NC}"
unset LD_LIBRARY_PATH
cd code
# 清理旧的构建产物，确保完全重新编译
echo "  - 清理旧构建..."
rm -f build/CMakeCache.txt
rm -f build/trace
# 重新配置和编译
echo "  - 配置CMake..."
if ! cmake -B build -DCMAKE_BUILD_TYPE=Release > "${PROJECT_ROOT}/${RESULT_DIR}/cmake.log" 2>&1; then
    echo -e "${RED}✗ CMake配置失败，查看 ${RESULT_DIR}/cmake.log${NC}"
    cat "${PROJECT_ROOT}/${RESULT_DIR}/cmake.log"
    exit 1
fi
echo "  - 编译..."
if cmake --build build --clean-first -j$(nproc) > "${PROJECT_ROOT}/${RESULT_DIR}/build.log" 2>&1; then
    # 验证可执行文件已生成且最新
    if [ ! -f "build/trace" ]; then
        echo -e "${RED}✗ 编译成功但未找到可执行文件 build/trace${NC}"
        exit 1
    fi
    TRACE_TIME=$(stat -c %Y build/trace)
    CURRENT_TIME=$(date +%s)
    TIME_DIFF=$((CURRENT_TIME - TRACE_TIME))
    if [ $TIME_DIFF -gt 60 ]; then
        echo -e "${RED}⚠ 警告: trace可执行文件不是最新的 (${TIME_DIFF}秒前)${NC}"
        echo -e "${RED}  可能存在缓存问题，请手动检查${NC}"
    fi
    echo -e "${GREEN}✓ 编译成功 (trace文件已更新)${NC}"
else
    echo -e "${RED}✗ 编译失败，查看 ${RESULT_DIR}/build.log${NC}"
    tail -20 "${PROJECT_ROOT}/${RESULT_DIR}/build.log"
    exit 1
fi
cd "$PROJECT_ROOT"
echo ""

# 2. 运行第三问
echo -e "${YELLOW}[2/6] 运行第三问...${NC}"
unset LD_LIBRARY_PATH
# 删除可能存在的旧输出文件
rm -f "${RESULT_DIR}/result.txt"
START_TIME=$(date +%s)
./code/build/trace \
    -layout dataset/trace_answer_checker/case/case1_small_layout.txt \
    -rule dataset/trace_answer_checker/Rule/public_small_rule3.txt \
    -output "${RESULT_DIR}/result.txt" \
    > "${RESULT_DIR}/run.log" 2>&1
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
# 验证输出文件已生成
if [ ! -f "${RESULT_DIR}/result.txt" ]; then
    echo -e "${RED}✗ 输出文件未生成${NC}"
    echo -e "${RED}  查看运行日志: ${RESULT_DIR}/run.log${NC}"
    tail -20 "${RESULT_DIR}/run.log"
    exit 1
fi
# 验证输出文件是最新的
RESULT_TIME=$(stat -c %Y "${RESULT_DIR}/result.txt")
TIME_DIFF=$((END_TIME - RESULT_TIME))
if [ $TIME_DIFF -gt 5 ]; then
    echo -e "${RED}⚠ 警告: 输出文件不是最新的 (${TIME_DIFF}秒前)${NC}"
fi
RESULT_SIZE=$(stat -c %s "${RESULT_DIR}/result.txt")
echo -e "${GREEN}✓ 运行完成 (耗时: ${RUNTIME}秒, 输出: ${RESULT_SIZE} 字节)${NC}"

# 提取关键信息
PHASE1_VISITED=$(grep "Phase 1:" "${RESULT_DIR}/run.log" | grep -oP 'visited \K[0-9]+' || echo "0")
PHASE1_HIGH=$(grep "Phase 1:" "${RESULT_DIR}/run.log" | grep -oP 'marked \K[0-9]+' || echo "0")
PHASE2_CUT=$(grep "Phase 2:" "${RESULT_DIR}/run.log" | grep -oP 'cut \K[0-9]+' || echo "0")
TOTAL_PIECES=$(grep "Total generated pieces:" "${RESULT_DIR}/run.log" | grep -oP ': \K[0-9]+' || echo "0")
VISITED_PIECES=$(grep "Total visited pieces:" "${RESULT_DIR}/run.log" | grep -oP ': \K[0-9]+' || echo "0")
UNVISITED_PIECES=$(grep "Unvisited pieces:" "${RESULT_DIR}/run.log" | grep -oP ': \K[0-9]+' || echo "0")
TOTAL_FOUND=$(grep "Found" "${RESULT_DIR}/run.log" | grep -oP 'Found \K[0-9]+' || echo "0")

echo ""
echo "程序输出统计:"
echo "  Phase 1 - 访问多边形: ${PHASE1_VISITED}"
echo "  Phase 1 - 高电平标记: ${PHASE1_HIGH}"
echo "  Phase 2 - AA切割数量: ${PHASE2_CUT}"
echo "  生成片段总数: ${TOTAL_PIECES}"
echo "  访问片段数量: ${VISITED_PIECES}"
echo "  未访问片段: ${UNVISITED_PIECES}"
echo "  总连通多边形: ${TOTAL_FOUND}"
echo ""

# 3. 统计各层多边形数量
echo -e "${YELLOW}[3/6] 统计各层多边形数量...${NC}"
cat > "${RESULT_DIR}/count_layers.py" << 'PYEOF'
import sys
from collections import defaultdict

layer_counts = defaultdict(int)
current_layer = None

with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('L'):
            current_layer = line
            layer_counts[current_layer] = 0
        elif current_layer and line.startswith('('):
            layer_counts[current_layer] += 1

for layer in sorted(layer_counts.keys()):
    print(f"{layer}: {layer_counts[layer]}")
PYEOF

python3 "${RESULT_DIR}/count_layers.py" "${RESULT_DIR}/result.txt" > "${RESULT_DIR}/layer_counts.txt"
echo -e "${GREEN}✓ 统计完成${NC}"
cat "${RESULT_DIR}/layer_counts.txt"
echo ""

# 4. 运行checker
echo -e "${YELLOW}[4/6] 运行checker测试...${NC}"
# 确保删除旧的临时文件
rm -f dataset/trace_answer_checker/test_q3_tmp.txt
# 复制最新结果
cp "${RESULT_DIR}/result.txt" dataset/trace_answer_checker/test_q3_tmp.txt
# 验证文件已复制且大小正确
if [ ! -f "dataset/trace_answer_checker/test_q3_tmp.txt" ]; then
    echo -e "${RED}✗ 结果文件复制失败${NC}"
    exit 1
fi
RESULT_SIZE=$(stat -c %s "${RESULT_DIR}/result.txt")
TMP_SIZE=$(stat -c %s "dataset/trace_answer_checker/test_q3_tmp.txt")
if [ "$RESULT_SIZE" != "$TMP_SIZE" ]; then
    echo -e "${RED}✗ 文件大小不匹配: 原始=${RESULT_SIZE}, 临时=${TMP_SIZE}${NC}"
    exit 1
fi
cd dataset/trace_answer_checker
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
./AnswerChecker -ans answer/public/small_txt/case1_small_q3_0912.txt -res test_q3_tmp.txt > "${PROJECT_ROOT}/${RESULT_DIR}/checker.log" 2>&1
cd "$PROJECT_ROOT"
echo -e "${GREEN}✓ Checker完成${NC}"
echo ""

# 5. 解析checker结果
echo -e "${YELLOW}[5/6] 解析checker结果...${NC}"
SCORE=$(grep "^score:" "${RESULT_DIR}/checker.log" | grep -oP ': \K[0-9.]+' || echo "0")

# 提取各层对比信息
cat > "${RESULT_DIR}/parse_checker.py" << 'PYEOF'
import sys
import re

with open(sys.argv[1], 'r') as f:
    content = f.read()

# 提取各层的poly数量
layer_pattern = r'\((\w+)\) compare finished'
size_pattern = r'poly1 size: (\d+), poly2 size: (\d+)'
status_pattern = r'(\w+): (pass|fail)'

layers = re.findall(layer_pattern, content)
sizes = re.findall(size_pattern, content)
statuses = re.findall(status_pattern, content)

print("层\t标准答案\t我们的结果\t差异\t\t差异率\t状态")
print("-" * 80)

for i, layer in enumerate(layers):
    if i < len(sizes):
        std = int(sizes[i][0])
        ours = int(sizes[i][1])
        diff = ours - std
        diff_pct = (diff / std * 100) if std > 0 else 0
        
        # 找到对应的状态
        status = "?"
        for s_layer, s_status in statuses:
            if s_layer == layer:
                status = s_status
                break
        
        status_icon = "✓" if status == "pass" else "✗"
        status_color = "pass" if status == "pass" else "fail"
        
        print(f"{layer}\t{std}\t\t{ours}\t\t{diff:+d}\t\t{diff_pct:+.1f}%\t{status_icon}")

print("")
print(f"总得分: {re.search(r'score: ([0-9.]+)', content).group(1) if re.search(r'score: ([0-9.]+)', content) else '0'}")

# 统计错误数量
errors = re.findall(r'compare error with (\d+) diffs', content)
if errors:
    total_errors = sum(int(e) for e in errors)
    print(f"总错误数: {total_errors}")
PYEOF

python3 "${RESULT_DIR}/parse_checker.py" "${RESULT_DIR}/checker.log" > "${RESULT_DIR}/layer_comparison.txt"
cat "${RESULT_DIR}/layer_comparison.txt"
echo ""

# 6. 生成详细报告
echo -e "${YELLOW}[6/6] 生成详细测试报告...${NC}"
cat > "${RESULT_DIR}/report.md" << MDEOF
# 第三问测试报告

**版本**: ${VERSION}  
**时间**: ${TIMESTAMP}  
**运行时间**: ${RUNTIME}秒  

## 1. 程序执行情况

### Phase 1 (高电平标记)
- 访问多边形数: ${PHASE1_VISITED}
- 高电平标记数: ${PHASE1_HIGH}

### Phase 2 (AA切割与扩展)
- AA切割数量: ${PHASE2_CUT}
- 生成片段总数: ${TOTAL_PIECES}
- 访问片段数量: ${VISITED_PIECES}
- 未访问片段数: ${UNVISITED_PIECES}

### 总结
- 总连通多边形: ${TOTAL_FOUND}

## 2. 各层多边形统计

\`\`\`
$(cat "${RESULT_DIR}/layer_counts.txt")
\`\`\`

## 3. Checker测试结果

**得分**: ${SCORE} / 60.00

### 各层详细对比

\`\`\`
$(cat "${RESULT_DIR}/layer_comparison.txt")
\`\`\`

## 4. 关键指标

| 指标 | 数值 |
|------|------|
| AA切割率 | ${PHASE2_CUT} / 45802 ($(python3 -c "print(f'{${PHASE2_CUT}/45802*100:.1f}%')")) |
| 片段访问率 | ${VISITED_PIECES} / ${TOTAL_PIECES} ($(python3 -c "print(f'{${VISITED_PIECES}/${TOTAL_PIECES}*100:.1f}%' if ${TOTAL_PIECES} > 0 else '0%')")) |
| L1层完成度 | $(grep "^L1:" "${RESULT_DIR}/layer_counts.txt" | awk '{printf "%.1f%%", ($2/113092)*100}') |

## 5. 问题诊断

$(python3 -c "
import sys
score = ${SCORE}
l1_count = $(grep "^L1:" "${RESULT_DIR}/layer_counts.txt" | cut -d' ' -f2)
aa_cut = ${PHASE2_CUT}

if score >= 55:
    print('✓ 得分优秀，算法表现良好')
elif score >= 40:
    print('△ 得分中等，需要进一步优化')
elif score >= 20:
    print('⚠ 得分较低，存在明显问题')
else:
    print('✗ 得分很低，需要重大改进')

print()

if l1_count < 100000:
    print(f'- L1层数量不足: {l1_count} < 113092')
    print('  建议: 检查AA切割逻辑和片段生成')
    
if aa_cut < 40000:
    print(f'- AA切割数量不足: {aa_cut} < 45802')
    print('  建议: 检查贯穿判据是否过于严格')
    
if ${UNVISITED_PIECES} > 0:
    print(f'- 存在未访问片段: {${UNVISITED_PIECES}}')
    print('  建议: 检查邻接关系构建和入口片段逻辑')
")

## 6. 文件位置

- 编译日志: \`${RESULT_DIR}/build.log\`
- 运行日志: \`${RESULT_DIR}/run.log\`
- 输出结果: \`${RESULT_DIR}/result.txt\`
- Checker日志: \`${RESULT_DIR}/checker.log\`
- 层统计: \`${RESULT_DIR}/layer_counts.txt\`
- 层对比: \`${RESULT_DIR}/layer_comparison.txt\`

---
生成时间: $(date)
MDEOF

echo -e "${GREEN}✓ 报告生成完成${NC}"
echo ""

# 显示总结
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  测试完成总结${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "版本: ${YELLOW}${VERSION}${NC}"
echo -e "得分: ${GREEN}${SCORE}${NC} / 60.00"
echo -e "AA切割: ${YELLOW}${PHASE2_CUT}${NC} / 45802"
echo -e "运行时间: ${YELLOW}${RUNTIME}${NC}秒"
echo ""
echo -e "详细报告: ${BLUE}${RESULT_DIR}/report.md${NC}"
echo ""

# 如果得分有提升，显示祝贺信息
if [ -f ".last_score" ]; then
    LAST_SCORE=$(cat .last_score)
    if (( $(echo "$SCORE > $LAST_SCORE" | bc -l) )); then
        IMPROVEMENT=$(echo "$SCORE - $LAST_SCORE" | bc)
        echo -e "${GREEN}🎉 恭喜！得分提升了 ${IMPROVEMENT} 分！${NC}"
        echo ""
    fi
fi
echo "$SCORE" > .last_score

# 清理临时文件
rm -f dataset/trace_answer_checker/test_q3_tmp.txt

echo -e "${GREEN}全部测试完成！${NC}"


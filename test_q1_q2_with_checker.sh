#!/bin/bash

set -e

PROJECT_ROOT="/root/autodl-tmp/EDA"
cd "$PROJECT_ROOT"

echo "========================================"
echo "  测试第一问和第二问（使用checker）"
echo "========================================"

# ============================================
# 测试第一问
# ============================================
echo ""
echo "【1】测试第一问"
echo "----------------------------------------"

# 每次都重新生成结果文件（使用最新编译的程序）
echo "生成第一问结果..."
rm -f test_q1_result.txt
unset LD_LIBRARY_PATH
./build/trace \
  -layout dataset/trace_answer_checker/case/case1_small_layout.txt \
  -rule dataset/trace_answer_checker/Rule/public_small_rule1.txt \
  -output test_q1_result.txt 2>&1 > /dev/null
echo "✓ 生成完成"

# 复制到checker目录
rm -f dataset/trace_answer_checker/test_q1_tmp.txt
cp test_q1_result.txt dataset/trace_answer_checker/test_q1_tmp.txt

# 验证文件已复制
if [ ! -f "dataset/trace_answer_checker/test_q1_tmp.txt" ]; then
    echo "✗ 结果文件复制失败"
    exit 1
fi

# 运行checker
cd dataset/trace_answer_checker
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
echo "运行checker..."
./AnswerChecker -ans answer/public/small_txt/case1_small_q1.txt -res test_q1_tmp.txt > "${PROJECT_ROOT}/test_q1_checker.log" 2>&1
cd "$PROJECT_ROOT"

# 解析结果
echo ""
echo "Checker结果:"
cat test_q1_checker.log | grep -E "score:|pass|fail|compare finished|poly1 size|poly2 size"

SCORE_Q1=$(grep "^score:" test_q1_checker.log | grep -oP ': \K[0-9.]+' || echo "0")
echo ""
echo "第一问得分: ${SCORE_Q1}"

# ============================================
# 测试第二问
# ============================================
echo ""
echo ""
echo "【2】测试第二问"
echo "----------------------------------------"

# 每次都重新生成结果文件（使用最新编译的程序）
echo "生成第二问结果..."
rm -f test_q2_result.txt
unset LD_LIBRARY_PATH
./build/trace \
  -layout dataset/trace_answer_checker/case/case1_small_layout.txt \
  -rule dataset/trace_answer_checker/Rule/public_small_rule2.txt \
  -output test_q2_result.txt 2>&1 > /dev/null
echo "✓ 生成完成"

# 复制到checker目录
rm -f dataset/trace_answer_checker/test_q2_tmp.txt
cp test_q2_result.txt dataset/trace_answer_checker/test_q2_tmp.txt

# 验证文件已复制
if [ ! -f "dataset/trace_answer_checker/test_q2_tmp.txt" ]; then
    echo "✗ 结果文件复制失败"
    exit 1
fi

# 运行checker
cd dataset/trace_answer_checker
export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH
echo "运行checker..."
./AnswerChecker -ans answer/public/small_txt/case1_small_q2.txt -res test_q2_tmp.txt > "${PROJECT_ROOT}/test_q2_checker.log" 2>&1
cd "$PROJECT_ROOT"

# 解析结果
echo ""
echo "Checker结果:"
cat test_q2_checker.log | grep -E "score:|pass|fail|compare finished|poly1 size|poly2 size"

SCORE_Q2=$(grep "^score:" test_q2_checker.log | grep -oP ': \K[0-9.]+' || echo "0")
echo ""
echo "第二问得分: ${SCORE_Q2}"

# ============================================
# 总结
# ============================================
echo ""
echo "========================================"
echo "  测试总结"
echo "========================================"
echo "第一问得分: ${SCORE_Q1}"
echo "第二问得分: ${SCORE_Q2}"
echo ""

if [ "$SCORE_Q1" == "60.00" ]; then
    echo "✅ 第一问：通过 (满分60.00)"
else
    echo "❌ 第一问：未通过（期望60.00，实际${SCORE_Q1}）"
fi

if [ "$SCORE_Q2" == "60.00" ]; then
    echo "✅ 第二问：通过 (满分60.00)"
else
    echo "❌ 第二问：未通过（期望60.00，实际${SCORE_Q2}）"
fi

echo ""
echo "详细日志:"
echo "  test_q1_checker.log"
echo "  test_q2_checker.log"

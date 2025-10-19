#pragma once
#include "types.hpp"
#include "gate.hpp"
#include <string>
#include <vector>

// 输出：按 Layer 分组，一行一个多边形；顶点逆时针
void write_result(const Layout& L, const std::vector<u64>& reachable, const GateCtx* gate, const std::string& out_path);

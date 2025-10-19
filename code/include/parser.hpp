#pragma once
#include "types.hpp"
#include <string>

Rule   parse_rule_file(const std::string& rule_path);
Layout parse_layout_file(const std::string& layout_path); // text -> memory

// NOTE: 可后续扩展为二进制缓存 + mmap 的版本：parse_layout_or_cache(...)

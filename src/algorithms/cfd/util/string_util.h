#pragma once

// see ../algorithms/cfd/LICENSE

#include <sstream>
#include <cstdarg>

[[maybe_unused]] std::string& LeftTrim(std::string &s);
[[maybe_unused]] std::string& RightTrim(std::string &s);
[[maybe_unused]] std::string& Trim(std::string &s);
[[maybe_unused]] std::string Concat(int count, ...);
[[maybe_unused]] std::string ConcatCsv(int count, ...);

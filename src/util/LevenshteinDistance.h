#pragma once

#include <string_view>

namespace util {

unsigned LevenshteinDistance(std::string_view l, std::string_view r);

}  // namespace util

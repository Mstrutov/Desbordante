#pragma once

#include <filesystem>
#include <optional>

#include "parser/csv_parser/csv_parser.h"

namespace util {

std::optional<char> ValidateSeparator(std::filesystem::path const& path, char separator);

}  // namespace util

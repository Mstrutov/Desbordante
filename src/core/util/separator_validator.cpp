#include "separator_validator.h"

#include <easylogging++.h>

namespace util {

std::optional<char> ValidateSeparator(std::filesystem::path const& path, char separator) {
    auto parser = std::make_unique<CSVParser>(path, separator, false);
    return parser->ValidateSeparator();
}

}  // namespace util

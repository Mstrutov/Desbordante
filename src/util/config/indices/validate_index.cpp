#include "util/config/indices/validate_index.h"

#include <stdexcept>

namespace util::config {

void ValidateIndex(IndexType value, size_t cols_count) {
    if (value >= cols_count) {
        throw std::invalid_argument(
                "Column index should be less than the number of columns in the dataset.");
    }
}

}  // namespace util::config

#pragma once

#include <cstddef>
#include <string>

#include "model/index.h"

namespace model::md {

struct ColumnMatch {
    Index left_col_index;
    Index right_col_index;
    std::string similarity_function_name;

    ColumnMatch(Index left_col_index, Index right_col_index,
                std::string similarity_function_name) noexcept
        : left_col_index(left_col_index),
          right_col_index(right_col_index),
          similarity_function_name(std::move(similarity_function_name)) {}
};

}  // namespace model::md

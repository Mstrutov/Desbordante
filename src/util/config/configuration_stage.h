#pragma once

#include <cstddef>

#include <enum.h>

namespace util::config {

// clang-format off
BETTER_ENUM(ConfigurationStage, size_t,
    load_data = 0,
    load_prepared_data,
    execute
)
// clang-format on

}  // namespace util::config

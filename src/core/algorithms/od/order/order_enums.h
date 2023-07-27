#pragma once

#include <enum.h>

namespace algos::order {

BETTER_ENUM(ValidityType, char,
    valid = 0,
    merge,
    swap
);

}    // namespace algos::order
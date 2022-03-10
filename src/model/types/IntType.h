#pragma once

#include "NumericType.h"

namespace model {

class IntType final : public NumericType<Int> {
public:
    IntType() noexcept : NumericType<Int>(TypeId::kInt) {}
};

}  // namespace model

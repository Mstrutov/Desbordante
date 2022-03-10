#pragma once

#include "NumericType.h"

namespace model {

class DoubleType final : public NumericType<Double> {
public:
    DoubleType() noexcept : NumericType<Double>(TypeId::kDouble) {}

    CompareResult Compare(std::byte const* l, std::byte const* r) const final {
        Double l_val = GetValue(l);
        Double r_val = GetValue(r);

        if (std::abs(l_val - r_val) < std::numeric_limits<Double>::epsilon()) {
            return CompareResult::kEqual;
        }

        if (l_val < r_val) {
            return CompareResult::kLess;
        }

        return CompareResult::kGreater;
    }
};

}  // namespace model

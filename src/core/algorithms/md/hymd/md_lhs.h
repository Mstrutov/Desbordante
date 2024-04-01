#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/utility/vector_double_hash.h"
#include "model/index.h"

namespace algos::hymd {
class MdLhs {
    std::vector<model::md::DecisionBoundary> values_;

public:
    MdLhs(std::size_t column_match_number) : values_(column_match_number, 0.0) {}

    std::vector<model::md::DecisionBoundary> const& GetValues() const noexcept {
        return values_;
    }

    model::md::DecisionBoundary& operator[](model::Index index) noexcept {
        return values_[index];
    }

    model::md::DecisionBoundary const& operator[](model::Index index) const noexcept {
        return values_[index];
    }

    std::size_t ColumnMatchNumber() const noexcept {
        return values_.size();
    }

    MdElement FindNextNonZero(model::Index start = 0) const noexcept {
        for (std::size_t const column_match_number = ColumnMatchNumber();
             start != column_match_number; ++start) {
            model::md::DecisionBoundary const bound = values_[start];
            if (bound != kLowestBound) return {start, bound};
        }
        return {start, kLowestBound};
    }

    bool IsNotEnd(MdElement element) const noexcept {
        return element.decision_boundary != kLowestBound;
    }

    bool IsEnd(MdElement element) const noexcept {
        return !IsNotEnd(element);
    }

    friend bool operator==(MdLhs const& lhs1, MdLhs const& lhs2) {
        return lhs1.GetValues() == lhs2.GetValues();
    }
};
}  // namespace algos::hymd

namespace std {
template <>
struct hash<algos::hymd::MdLhs> {
    std::size_t operator()(algos::hymd::MdLhs const& p) const noexcept {
        return hash<algos::hymd::DecisionBoundaryVector>{}(p.GetValues());
    }
};
}  // namespace std

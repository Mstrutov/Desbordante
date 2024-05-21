#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/utility/vector_double_hash.h"
#include "model/index.h"
#include "util/py_tuple_hash.h"

namespace algos::hymd {
struct LhsNode {
    model::Index child_array_index;
    model::md::DecisionBoundary decision_boundary;

    friend bool operator==(LhsNode const& l, LhsNode const& r) {
        return l.child_array_index == r.child_array_index &&
               l.decision_boundary == r.decision_boundary;
    }
};

class MdLhs {
    using Nodes = std::vector<LhsNode>;
    Nodes values_;

public:
    using iterator = Nodes::const_iterator;

    MdLhs(std::size_t column_match_number) {
        values_.reserve(column_match_number);
    }

    model::md::DecisionBoundary& AddNext(model::Index child_array_index) {
        return values_.emplace_back(child_array_index).decision_boundary;
    }

    void RemoveLast() {
        values_.pop_back();
    }

    iterator begin() const noexcept {
        return values_.begin();
    }

    iterator end() const noexcept {
        return values_.end();
    }

    friend bool operator==(MdLhs const& lhs1, MdLhs const& lhs2) {
        return lhs1.values_ == lhs2.values_;
    }

    std::size_t Cardinality() const noexcept {
        return values_.size();
    }

    bool IsEmpty() const noexcept {
        return values_.empty();
    }
};

}  // namespace algos::hymd

namespace std {
template <>
struct hash<algos::hymd::MdLhs> {
    std::size_t operator()(algos::hymd::MdLhs const& p) const noexcept {
        using model::Index, model::md::DecisionBoundary;
        util::PyTupleHash main_hasher(p.Cardinality());
        for (auto const& [child_array_index, bound] : p) {
            util::PyTupleHash pair_hasher(2);
            pair_hasher.AppendHash(std::hash<Index>{}(child_array_index));
            pair_hasher.AppendHash(std::hash<DecisionBoundary>{}(bound));
            main_hasher.AppendHash(pair_hasher.GetResult());
        }
        return main_hasher.GetResult();
    }
};
}  // namespace std

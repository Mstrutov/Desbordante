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
    std::size_t column_match_number_;

public:
    using iterator = Nodes::const_iterator;

    MdLhs(std::size_t column_match_number) : column_match_number_(column_match_number) {
        values_.reserve(column_match_number_);
    }

    model::md::DecisionBoundary& AddNext(model::Index child_array_index) {
        return values_.emplace_back(child_array_index).decision_boundary;
    }

    void RemoveLast() {
        values_.pop_back();
    }

    std::vector<model::md::DecisionBoundary> ToBoundVec() const {
        std::vector<model::md::DecisionBoundary> bound_vec(ColumnMatchNumber(), kLowestBound);
        model::Index cur_index = 0;
        for (auto const& [child_array_index, bound] : values_) {
            cur_index += child_array_index;
            bound_vec[cur_index] = bound;
            ++cur_index;
        }
        return bound_vec;
    }

    model::md::DecisionBoundary operator[](model::Index index) const {
        return ToBoundVec()[index];
    }

    std::size_t ColumnMatchNumber() const noexcept {
        return column_match_number_;
    }

    MdElement FindNextNonZero(model::Index start = 0) const noexcept {
        auto bounds = ToBoundVec();
        for (std::size_t const column_match_number = ColumnMatchNumber();
             start != column_match_number; ++start) {
            model::md::DecisionBoundary const bound = bounds[start];
            if (bound != kLowestBound) return {start, bound};
        }
        return {start, kLowestBound};
    }

    iterator begin() const noexcept {
        return values_.begin();
    }

    iterator end() const noexcept {
        return values_.end();
    }

    model::Index ToIndex(iterator iter) const noexcept {
        if (iter == end()) return ColumnMatchNumber();
        model::Index index = 0;
        auto cur_iter = begin();
        for (; cur_iter != iter; ++cur_iter, ++index) {
            index += cur_iter->child_array_index;
        }
        index += cur_iter->child_array_index;
        return index;
    }

    iterator FindIter(model::Index index) const noexcept {
        iterator iter = begin(), end_iter = end();
        for (model::Index cur_index = 0; iter != end_iter; ++iter, ++cur_index) {
            cur_index += iter->child_array_index;
            if (cur_index >= index) return iter;
        }
        return iter;
    }

    bool IsNotEnd(MdElement element) const noexcept {
        return element.decision_boundary != kLowestBound;
    }

    bool IsEnd(MdElement element) const noexcept {
        return !IsNotEnd(element);
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

    model::Index GetColumnMatchIndex(iterator iter, model::Index child_index) const noexcept {
        model::Index index = 0;
        for (auto it = begin(); it != iter; ++it) {
            index += it->child_array_index;
            ++index;
        }
        return index + child_index;
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

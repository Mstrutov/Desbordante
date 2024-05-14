#include "algorithms/md/hymd/lattice/cardinality/one_by_one_min_picker.h"

#include <algorithm>
#include <cassert>

namespace algos::hymd::lattice::cardinality {

void OneByOnePicker::NewBatch(std::size_t elements) {
    currently_picked_.clear();
    currently_picked_.reserve(elements);
}

// Assuming [gen_it, gen_end) does not specialize [spec_it, spec_end), is the former a
// generalization of the latter?
bool OneByOnePicker::IsGeneralization(MdLhs::iterator gen_it, MdLhs::iterator spec_it,
                                      MdLhs::iterator const gen_end,
                                      MdLhs::iterator const spec_end) {
    using model::Index;
    auto inc_until_eq = [](MdLhs::iterator& gen_it, MdLhs::iterator& spec_it,
                           MdLhs::iterator const spec_end) {
        Index cur_spec_index = spec_it->child_array_index;
        Index const gen_index = gen_it->child_array_index;
        // Incomparable, spec has no classifier for the considered gen's column match.
        if (cur_spec_index > gen_index) return true;
        if (cur_spec_index == gen_index) {
            // Incomparable, spec's classifier matches more record pairs.
            if (gen_it->decision_boundary > spec_it->decision_boundary) return true;
            ++gen_it;
            ++spec_it;
            return false;
        }
        while (true) {
            ++spec_it;
            // Incomparable, generalization's node is after the last specialization's node.
            if (spec_it == spec_end) return true;
            ++cur_spec_index;
            cur_spec_index += spec_it->child_array_index;
            // Incomparable, spec has no classifier for the considered gen's column match.
            if (cur_spec_index > gen_index) return true;
            if (cur_spec_index == gen_index) {
                // Incomparable, spec's classifier matches more record pairs.
                if (gen_it->decision_boundary > spec_it->decision_boundary) return true;
                ++gen_it;
                ++spec_it;
                return false;
            }
        }
    };
    while (true) {
        bool is_incomparable = inc_until_eq(gen_it, spec_it, spec_end);
        if (is_incomparable) return false;
        // Fewer classifiers in gen, gen generalizes spec.
        if (gen_it == gen_end) return true;
        // More classifiers in gen, gen is incomparable with spec.
        if (spec_it == spec_end) return false;
    }
}

auto OneByOnePicker::CompareLhss(MdLhs const& cur, MdLhs const& prev) -> ComparisonResult {
    auto cur_it = cur.begin();
    auto const cur_end = cur.end();
    auto prev_it = prev.begin();
    auto const prev_end = prev.end();
    while (cur_it != cur_end) {     // cur still has LHS elements.
        if (prev_it == prev_end) {  // prev has no more LHS elements, prev generalizes cur.
            return ComparisonResult::Specialization;
        }
        auto const& [cur_child_index, cur_bound] = *cur_it;
        auto const& [prev_child_index, prev_bound] = *prev_it;
        if (cur_child_index > prev_child_index ||
            (cur_child_index == prev_child_index &&
             cur_bound < prev_bound)) {  // prev cannot be a generalization of cur.
            if (IsGeneralization(cur_it, prev_it, cur_end, prev_end))
                return ComparisonResult::Generalization;
            return ComparisonResult::Incomparable;
        }
        if (cur_child_index < prev_child_index ||
            (cur_child_index == prev_child_index &&
             prev_bound < cur_bound)) {  // cur cannot be a generalization of prev.
            if (IsGeneralization(prev_it, cur_it, prev_end, cur_end)) {
                return ComparisonResult::Specialization;
            }
            return ComparisonResult::Incomparable;
        }
        ++cur_it;
        ++prev_it;
    }
    // Assuming all LHSs are not equal.
    assert(prev_it != prev_end);
    // cur has no more LHS elements, cur generalizes prev.
    return ComparisonResult::Generalization;
}

void OneByOnePicker::AddGeneralizations(MdLattice::MdVerificationMessenger& messenger,
                                        boost::dynamic_bitset<>& considered_indices) {
    MdLhs const& lhs_cur = messenger.GetLhs();
    assert(!lhs_cur.IsEmpty() || currently_picked_.empty());
    for (ValidationInfo& prev_info : currently_picked_) {
        MdLhs const& lhs_prev = prev_info.messenger->GetLhs();
        boost::dynamic_bitset<>& indices_prev = prev_info.rhs_indices;
        switch (CompareLhss(lhs_cur, lhs_prev)) {
            case ComparisonResult::Specialization:
                considered_indices -= indices_prev;
                if (considered_indices.none()) return;
                break;
            case ComparisonResult::Generalization:
                indices_prev -= considered_indices;
                break;
            case ComparisonResult::Incomparable:
                break;
            default:
                assert(false);
                __builtin_unreachable();
        }
    }
    assert(!considered_indices.none());
    currently_picked_.emplace_back(&messenger, std::move(considered_indices));
}

std::vector<ValidationInfo> OneByOnePicker::GetAll() noexcept {
    static_assert(kNeedsEmptyRemoval,
                  "This picker needs post-processing to remove candidates with empty RHS indices");
    return std::move(currently_picked_);
}

}  // namespace algos::hymd::lattice::cardinality

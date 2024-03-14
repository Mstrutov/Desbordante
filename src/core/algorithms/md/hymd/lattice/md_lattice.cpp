#include "algorithms/md/hymd/lattice/md_lattice.h"

#include <algorithm>
#include <cassert>

#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/utility/get_first_non_zero_index.h"
#include "algorithms/md/hymd/utility/set_for_scope.h"
#include "util/erase_if_replace.h"

namespace {
using namespace algos::hymd::lattice;
using namespace algos::hymd;
using model::md::DecisionBoundary, model::Index;
using utility::GetFirstNonZeroIndex;

bool NotEmpty(DecisionBoundaryVector const& rhs_bounds) {
    auto not_lowest = [](DecisionBoundary bound) { return bound != kLowestBound; };
    return std::any_of(rhs_bounds.begin(), rhs_bounds.end(), not_lowest);
}
}  // namespace

namespace algos::hymd::lattice {

// TODO: remove recursion
MdLattice::MdLattice(std::size_t column_matches_size, SingleLevelFunc single_level_func,
                     std::vector<ColumnMatchInfo> const& column_matches_info,
                     bool prune_nondisjoint)
    : column_matches_size_(column_matches_size),
      md_root_(DecisionBoundaryVector(column_matches_size_, 1.0)),
      support_root_(column_matches_size_),
      get_single_level_(std::move(single_level_func)),
      column_matches_info_(&column_matches_info),
      prune_nondisjoint_(prune_nondisjoint) {}

std::optional<DecisionBoundary> MdLattice::SpecializeOneLhs(Index col_match_index,
                                                            DecisionBoundary lhs_bound) const {
    std::vector<DecisionBoundary> const& decision_bounds =
            (*column_matches_info_)[col_match_index].similarity_info.lhs_bounds;
    auto end_bounds = decision_bounds.end();
    auto upper = std::upper_bound(decision_bounds.begin(), end_bounds, lhs_bound);
    if (upper == end_bounds) return std::nullopt;
    return *upper;
}

void MdLattice::Specialize(DecisionBoundaryVector& lhs_bounds,
                           DecisionBoundaryVector const& specialize_past, Rhss const& rhss) {
    auto specialize_all_lhs = [this, &lhs_bounds, it_begin = rhss.begin(), it_end = rhss.end(),
                               &specialize_past](auto handle_same_lhs_as_rhs) {
        for (Index lhs_spec_index = 0; lhs_spec_index < column_matches_size_; ++lhs_spec_index) {
            std::optional<DecisionBoundary> const specialized_lhs_bound =
                    SpecializeOneLhs(lhs_spec_index, specialize_past[lhs_spec_index]);
            if (!specialized_lhs_bound.has_value()) continue;
            auto context = utility::SetForScope(lhs_bounds[lhs_spec_index], *specialized_lhs_bound);
            if (IsUnsupported(lhs_bounds)) continue;

            for (auto it = it_begin; it != it_end; ++it) {
                if (it->index == lhs_spec_index) {
                    handle_same_lhs_as_rhs(*it, *specialized_lhs_bound);
                    for (++it; it != it_end; ++it) {
                        AddIfMinimal(lhs_bounds, *it);
                    }
                    break;
                }
                AddIfMinimal(lhs_bounds, *it);
            }
        }
    };
    if (prune_nondisjoint_) {
        specialize_all_lhs([](...) {});
    } else {
        specialize_all_lhs(
                [this, &lhs_bounds](MdElement old_rhs, DecisionBoundary specialized_lhs_bound) {
                    if (old_rhs.decision_boundary > specialized_lhs_bound) {
                        AddIfMinimal(lhs_bounds, old_rhs);
                    }
                });
    }
}

void MdLattice::MdRefiner::Refine() {
    for (auto upd_iter = invalidated_.UpdateIterBegin(), upd_end = invalidated_.UpdateIterEnd();
         upd_iter != upd_end; ++upd_iter) {
        auto const& [rhs_index, new_bound] = *upd_iter;
        DecisionBoundary& md_rhs_bound_ref = (*rhs_)[rhs_index];
        md_rhs_bound_ref = kLowestBound;
        // trivial
        if (new_bound <= lhs_[rhs_index]) continue;
        // not minimal
        if (lattice_->HasGeneralization(lhs_, *upd_iter)) continue;
        md_rhs_bound_ref = new_bound;
    }
    lattice_->Specialize(lhs_, *sim_, invalidated_.GetInvalidated());
}

void MdLattice::TryAddRefiner(std::vector<MdRefiner>& found, DecisionBoundaryVector& rhs,
                              SimilarityVector const& similarity_vector,
                              DecisionBoundaryVector const& cur_node_lhs_bounds) {
    utility::InvalidatedRhss invalidated;
    for (Index i = 0; i != column_matches_size_; ++i) {
        DecisionBoundary sim_bound = similarity_vector[i];
        DecisionBoundary rhs_bound = rhs[i];
        if (sim_bound >= rhs_bound) continue;
        invalidated.PushBack({i, rhs_bound}, sim_bound);
        for (++i; i != column_matches_size_; ++i) {
            DecisionBoundary sim_bound = similarity_vector[i];
            DecisionBoundary rhs_bound = rhs[i];
            if (sim_bound >= rhs_bound) continue;
            invalidated.PushBack({i, rhs_bound}, sim_bound);
        }
        found.emplace_back(this, &similarity_vector, cur_node_lhs_bounds, &rhs,
                           std::move(invalidated));
        break;
    }
}

void MdLattice::CollectRefinersForViolated(MdNode& cur_node, std::vector<MdRefiner>& found,
                                           DecisionBoundaryVector& cur_node_lhs_bounds,
                                           SimilarityVector const& similarity_vector,
                                           Index cur_node_index) {
    TryAddRefiner(found, cur_node.rhs_bounds, similarity_vector, cur_node_lhs_bounds);

    MdNodeChildren& children = cur_node.children;
    std::size_t const child_array_size = children.size();
    for (Index child_array_index = FindFirstNonEmptyIndex(children, 0);
         child_array_index != child_array_size;
         child_array_index = FindFirstNonEmptyIndex(children, child_array_index + 1)) {
        Index const next_node_index = cur_node_index + child_array_index;
        DecisionBoundary& cur_lhs_bound = cur_node_lhs_bounds[next_node_index];
        DecisionBoundary const sim_vec_sim = similarity_vector[next_node_index];
        for (auto& [generalization_bound, node] : *children[child_array_index]) {
            if (generalization_bound > sim_vec_sim) break;
            cur_lhs_bound = generalization_bound;
            CollectRefinersForViolated(node, found, cur_node_lhs_bounds, similarity_vector,
                                       next_node_index + 1);
        }
        cur_lhs_bound = kLowestBound;
    }
}

auto MdLattice::CollectRefinersForViolated(SimilarityVector const& similarity_vector)
        -> std::vector<MdRefiner> {
    std::vector<MdRefiner> found;
    DecisionBoundaryVector current_lhs(column_matches_size_, kLowestBound);
    CollectRefinersForViolated(md_root_, found, current_lhs, similarity_vector, 0);
    return found;
}

void MdLattice::MdVerificationMessenger::MarkUnsupported() {
    // TODO: specializations can be removed from the MD lattice. If not worth it, removing just
    // this node and its children should be cheap. Though, destructors also take time.
    lattice_->MarkUnsupported(lhs_);
}

void MdLattice::MdVerificationMessenger::LowerAndSpecialize(
        utility::InvalidatedRhss const& invalidated) {
    for (auto upd_iter = invalidated.UpdateIterBegin(), upd_end = invalidated.UpdateIterEnd();
         upd_iter != upd_end; ++upd_iter) {
        auto const& [rhs_index, new_bound] = *upd_iter;
        (*rhs_)[rhs_index] = new_bound;
    }
    lattice_->Specialize(lhs_, lhs_, invalidated.GetInvalidated());
}

void MdLattice::AddNewMinimal(MdNode& cur_node, DecisionBoundaryVector const& lhs_bounds,
                              MdElement rhs, Index cur_node_index) {
    assert(IsEmpty(cur_node.children));
    assert(!NotEmpty(cur_node.rhs_bounds));
    MdNode* cur_node_ptr = &cur_node;
    for (Index next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != column_matches_size_; cur_node_index = next_node_index + 1,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        std::size_t const child_array_index = next_node_index - cur_node_index;
        std::size_t const next_child_array_size = column_matches_size_ - next_node_index;
        cur_node_ptr = &cur_node_ptr->children[child_array_index]
                                .emplace()
                                .try_emplace(lhs_bounds[next_node_index], column_matches_size_,
                                             next_child_array_size)
                                .first->second;
    }
    cur_node_ptr->rhs_bounds[rhs.index] = rhs.decision_boundary;
    UpdateMaxLevel(lhs_bounds);
}

bool MdLattice::HasLhsGeneralization(MdNode const& node, DecisionBoundaryVector const& lhs_bounds,
                                     MdElement rhs, Index const node_index,
                                     Index const start_index) const {
    for (Index next_node_index = GetFirstNonZeroIndex(lhs_bounds, start_index);
         next_node_index != column_matches_size_;
         next_node_index = GetFirstNonZeroIndex(lhs_bounds, next_node_index + 1)) {
        Index const child_array_index = next_node_index - node_index;
        MdOptionalChild const& optional_child = node.children[child_array_index];
        if (!optional_child.has_value()) continue;
        DecisionBoundary const generalization_bound_limit = lhs_bounds[next_node_index];
        for (auto const& [generalization_bound, node] : *optional_child) {
            if (generalization_bound > generalization_bound_limit) break;
            if (HasGeneralization(node, lhs_bounds, rhs, next_node_index + 1)) return true;
        }
    }
    return false;
}

void MdLattice::UpdateMaxLevel(DecisionBoundaryVector const& lhs_bounds) {
    std::size_t level = 0;
    for (Index i = 0; i != column_matches_size_; ++i) {
        DecisionBoundary cur_bound = lhs_bounds[i];
        if (cur_bound == kLowestBound) continue;
        level += get_single_level_(cur_bound, i);
    }
    if (level > max_level_) max_level_ = level;
}

class MdLattice::GeneralizationChecker {
private:
    MdNode* node_ = nullptr;
    DecisionBoundary* cur_node_rhs_ptr_ = nullptr;
    MdElement const rhs_;

public:
    GeneralizationChecker(MdElement rhs) noexcept : rhs_(rhs) {}

    bool SetAndCheck(MdNode& node) noexcept {
        node_ = &node;
        cur_node_rhs_ptr_ = &node_->rhs_bounds[rhs_.index];
        return *cur_node_rhs_ptr_ >= rhs_.decision_boundary;
    }

    MdNode& CurNode() noexcept {
        return *node_;
    }

    MdNodeChildren& Children() noexcept {
        return node_->children;
    }

    DecisionBoundary IncomingBound() const noexcept {
        return rhs_.decision_boundary;
    }

    Index IncomingIndex() const noexcept {
        return rhs_.index;
    }

    void SetBoundOnCurrent() noexcept {
        *cur_node_rhs_ptr_ = rhs_.decision_boundary;
    }

    MdElement GetIncomingRhs() const noexcept {
        return rhs_;
    }
};

auto MdLattice::ReturnNextNode(DecisionBoundaryVector const& lhs_bounds,
                               GeneralizationChecker& checker, Index cur_node_index,
                               Index next_node_index) -> MdNode* {
    MdElement const rhs = checker.GetIncomingRhs();
    Index const child_array_index = next_node_index - cur_node_index;
    auto [boundary_mapping, is_first_arr] = TryEmplaceChild(checker.Children(), child_array_index);
    DecisionBoundary const next_lhs_bound = lhs_bounds[next_node_index];
    std::size_t const next_child_array_size = column_matches_size_ - next_node_index;
    if (is_first_arr) [[unlikely]] {
        MdNode& new_node =
                boundary_mapping
                        .try_emplace(next_lhs_bound, column_matches_size_, next_child_array_size)
                        .first->second;
        AddNewMinimal(new_node, lhs_bounds, rhs, next_node_index + 1);
        return nullptr;
    }

    auto it = boundary_mapping.begin();
    for (auto end = boundary_mapping.end(); it != end; ++it) {
        auto& [generalization_bound, node] = *it;
        if (generalization_bound > next_lhs_bound) break;
        if (generalization_bound == next_lhs_bound) return &it->second;
        if (HasGeneralization(node, lhs_bounds, rhs, next_node_index + 1)) return nullptr;
    }
    using std::forward_as_tuple;
    MdNode& new_node =
            boundary_mapping
                    .emplace_hint(it, std::piecewise_construct, forward_as_tuple(next_lhs_bound),
                                  forward_as_tuple(column_matches_size_, next_child_array_size))
                    ->second;
    AddNewMinimal(new_node, lhs_bounds, rhs, next_node_index + 1);
    return nullptr;
}

void MdLattice::AddIfMinimal(DecisionBoundaryVector const& lhs_bounds, MdElement rhs) {
    // TODO: use info about where the LHS was specialized from to reduce generalization checks.
    // When an MD is inferred from, it is not a specialization of any other MD in the lattice, so
    // its LHS is not a specialization of any other. The only LHSs we have to check after a
    // specialization are those that are generalizations of the new LHS but not the old one.
    assert(!IsUnsupported(lhs_bounds));

    auto checker = GeneralizationChecker(rhs);
    if (checker.SetAndCheck(md_root_)) return;

    for (Index cur_node_index = 0,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != column_matches_size_; cur_node_index = next_node_index + 1,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        if (HasLhsGeneralization(checker.CurNode(), lhs_bounds, rhs, cur_node_index,
                                 next_node_index + 1))
            return;

        MdNode* next_node = ReturnNextNode(lhs_bounds, checker, cur_node_index, next_node_index);
        if (!next_node) return;
        if (checker.SetAndCheck(*next_node)) return;
    }
    // Correct. This method is only used when inferring, so we must get a maximum, which is
    // enforced with either of `if (cur_node_rhs_bound >= rhs_bound) return false;` above, no need
    // for an additional check.
    // I believe Metanome implemented this method incorrectly (they used Math.min instead of
    // Math.max). However, before validating an MD, they set the rhs bound to 1.0, so that error has
    // no effect on the final result. If this boundary is set correctly (like here), we don't have
    // to set the RHS bound to 1.0 before validating and can start with the value that was
    // originally in the lattice, potentially yielding more recommendations, since the only way that
    // value could be not equal to 1.0 at the start of validation is if it was lowered by some
    // record pair.
    checker.SetBoundOnCurrent();
}

bool MdLattice::HasGeneralization(MdNode const& node, DecisionBoundaryVector const& lhs_bounds,
                                  MdElement rhs, Index const cur_node_index) const {
    if (node.rhs_bounds[rhs.index] >= rhs.decision_boundary) return true;
    if (HasLhsGeneralization(node, lhs_bounds, rhs, cur_node_index, cur_node_index)) return true;
    return false;
}

bool MdLattice::HasGeneralization(DecisionBoundaryVector const& lhs_bounds, MdElement rhs) const {
    return HasGeneralization(md_root_, lhs_bounds, rhs, 0);
}

void MdLattice::RaiseInterestingnessBounds(
        MdNode const& cur_node, DecisionBoundaryVector const& lhs_bounds,
        std::vector<DecisionBoundary>& cur_interestingness_bounds, Index const this_node_index,
        std::vector<Index> const& indices) const {
    {
        std::size_t const indices_size = indices.size();
        for (Index i = 0; i < indices_size; ++i) {
            DecisionBoundary const this_node_bound = cur_node.rhs_bounds[indices[i]];
            DecisionBoundary& cur = cur_interestingness_bounds[i];
            if (this_node_bound > cur) {
                cur = this_node_bound;
                // The original paper mentions checking for the case where all decision bounds are
                // 1.0, but if such a situation occurs for any one RHS, and the generalization with
                // that RHS happens to be valid on the data, it would make inference from record
                // pairs give an incorrect result, meaning the algorithm is incorrect.
                // However, it is possible to stop decreasing when the bound's index in the list of
                // natural decision boundaries is exactly one less than the RHS bound's index.
                // TODO: abort traversal as above.
                assert(this_node_bound != 1.0);
            }
        }
    }

    std::size_t const col_match_number = lhs_bounds.size();
    for (Index next_node_index = GetFirstNonZeroIndex(lhs_bounds, this_node_index);
         next_node_index != col_match_number;
         next_node_index = GetFirstNonZeroIndex(lhs_bounds, next_node_index + 1)) {
        Index const child_array_index = next_node_index - this_node_index;
        MdOptionalChild const& optional_child = cur_node.children[child_array_index];
        if (!optional_child.has_value()) continue;
        DecisionBoundary const generalization_bound_limit = lhs_bounds[next_node_index];
        for (auto const& [generalization_bound, node] : *optional_child) {
            if (generalization_bound > generalization_bound_limit) break;
            RaiseInterestingnessBounds(node, lhs_bounds, cur_interestingness_bounds,
                                       next_node_index + 1, indices);
        }
    }
}

std::vector<DecisionBoundary> MdLattice::GetRhsInterestingnessBounds(
        DecisionBoundaryVector const& lhs_bounds, std::vector<Index> const& indices) const {
    std::vector<DecisionBoundary> interestingness_bounds;
    interestingness_bounds.reserve(indices.size());
    for (Index index : indices) {
        DecisionBoundary const lhs_bound = lhs_bounds[index];
        assert(lhs_bound != 1.0);
        interestingness_bounds.push_back(lhs_bound);
    }
    RaiseInterestingnessBounds(md_root_, lhs_bounds, interestingness_bounds, 0, indices);
    return interestingness_bounds;
}

void MdLattice::GetLevel(MdNode& cur_node, std::vector<MdVerificationMessenger>& collected,
                         DecisionBoundaryVector& cur_node_lhs_bounds, Index const cur_node_index,
                         std::size_t const level_left) {
    DecisionBoundaryVector& rhs_bounds = cur_node.rhs_bounds;
    MdNodeChildren& children = cur_node.children;
    if (level_left == 0) {
        if (NotEmpty(rhs_bounds)) collected.emplace_back(this, cur_node_lhs_bounds, &rhs_bounds);
        return;
    }
    std::size_t const child_array_size = children.size();
    for (Index child_array_index = FindFirstNonEmptyIndex(children, 0);
         child_array_index != child_array_size;
         child_array_index = FindFirstNonEmptyIndex(children, child_array_index + 1)) {
        Index const next_node_index = cur_node_index + child_array_index;
        DecisionBoundary& next_lhs_bound = cur_node_lhs_bounds[next_node_index];
        for (auto& [boundary, node] : *children[child_array_index]) {
            assert(boundary > kLowestBound);
            std::size_t const single = get_single_level_(next_node_index, boundary);
            if (single > level_left) break;
            next_lhs_bound = boundary;
            GetLevel(node, collected, cur_node_lhs_bounds, next_node_index + 1,
                     level_left - single);
        }
        next_lhs_bound = kLowestBound;
    }
}

auto MdLattice::GetLevel(std::size_t const level) -> std::vector<MdVerificationMessenger> {
    // TODO: traverse both simultaneously.
    std::vector<MdVerificationMessenger> collected;
    DecisionBoundaryVector current_lhs(column_matches_size_, kLowestBound);
    GetLevel(md_root_, collected, current_lhs, 0, level);
    util::EraseIfReplace(collected, [this](MdVerificationMessenger const& messenger) {
        return IsUnsupported(messenger.GetLhs());
    });
    return collected;
}

void MdLattice::GetAll(MdNode& cur_node, std::vector<MdLatticeNodeInfo>& collected,
                       DecisionBoundaryVector& this_node_lhs_bounds, Index const this_node_index) {
    MdNodeChildren& children = cur_node.children;
    DecisionBoundaryVector& rhs_bounds = cur_node.rhs_bounds;
    if (NotEmpty(rhs_bounds)) collected.emplace_back(this_node_lhs_bounds, &rhs_bounds);
    std::size_t const child_array_size = children.size();
    for (Index child_array_index = FindFirstNonEmptyIndex(children, 0);
         child_array_index != child_array_size;
         child_array_index = FindFirstNonEmptyIndex(children, child_array_index + 1)) {
        Index const next_node_index = this_node_index + child_array_index;
        DecisionBoundary& next_lhs_bound = this_node_lhs_bounds[next_node_index];
        for (auto& [boundary, node] : *children[child_array_index]) {
            next_lhs_bound = boundary;
            GetAll(node, collected, this_node_lhs_bounds, next_node_index + 1);
        }
        next_lhs_bound = kLowestBound;
    }
}

std::vector<MdLatticeNodeInfo> MdLattice::GetAll() {
    std::vector<MdLatticeNodeInfo> collected;
    DecisionBoundaryVector current_lhs(column_matches_size_, kLowestBound);
    GetAll(md_root_, collected, current_lhs, 0);
    util::EraseIfReplace(collected, [this](MdLatticeNodeInfo const& node_info) {
        return IsUnsupported(node_info.lhs_bounds);
    });
    return collected;
}

bool MdLattice::IsUnsupported(SupportNode const& cur_node, DecisionBoundaryVector const& lhs_bounds,
                              Index this_node_index) const {
    if (cur_node.is_unsupported) return true;
    SupportNodeChildren const& children = cur_node.children;
    std::size_t const child_array_size = children.size();
    for (Index child_array_index = FindFirstNonEmptyIndex(children, 0);
         child_array_index != child_array_size;
         child_array_index = FindFirstNonEmptyIndex(children, child_array_index + 1)) {
        Index const next_node_index = this_node_index + child_array_index;
        DecisionBoundary const generalization_boundary_limit = lhs_bounds[next_node_index];
        for (auto const& [generalization_boundary, node] : *children[child_array_index]) {
            if (generalization_boundary > generalization_boundary_limit) break;
            if (IsUnsupported(node, lhs_bounds, next_node_index + 1)) return true;
        }
    }
    return false;
}

void MdLattice::MarkNewLhs(SupportNode& cur_node, DecisionBoundaryVector const& lhs_bounds,
                           Index cur_node_index) {
    assert(IsEmpty(cur_node.children));
    std::size_t const col_match_number = lhs_bounds.size();
    SupportNode* cur_node_ptr = &cur_node;
    for (Index next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != col_match_number; cur_node_index = next_node_index + 1,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        std::size_t const child_array_index = next_node_index - cur_node_index;
        std::size_t const next_child_array_size = col_match_number - next_node_index;
        cur_node_ptr = &cur_node_ptr->children[child_array_index]
                                .emplace()
                                .try_emplace(lhs_bounds[next_node_index], next_child_array_size)
                                .first->second;
    }
    cur_node_ptr->is_unsupported = true;
}

bool MdLattice::IsUnsupported(DecisionBoundaryVector const& lhs_bounds) const {
    return IsUnsupported(support_root_, lhs_bounds, 0);
}

void MdLattice::MarkUnsupported(DecisionBoundaryVector const& lhs_bounds) {
    std::size_t const col_match_number = lhs_bounds.size();
    SupportNode* cur_node_ptr = &support_root_;
    for (Index cur_node_index = 0,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != col_match_number; cur_node_index = next_node_index + 1,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        DecisionBoundary const next_bound = lhs_bounds[next_node_index];
        Index const child_array_index = next_node_index - cur_node_index;
        std::size_t const next_child_array_size = col_match_number - next_node_index;
        auto [boundary_map, is_first_arr] =
                TryEmplaceChild(cur_node_ptr->children, child_array_index);
        if (is_first_arr) {
            SupportNode& new_node =
                    boundary_map.try_emplace(next_bound, next_child_array_size).first->second;
            MarkNewLhs(new_node, lhs_bounds, next_node_index + 1);
            return;
        }
        auto [it_map, is_first_map] = boundary_map.try_emplace(next_bound, next_child_array_size);
        SupportNode& next_node = it_map->second;
        if (is_first_map) {
            MarkNewLhs(next_node, lhs_bounds, next_node_index + 1);
            return;
        }
        cur_node_ptr = &next_node;
    }
    // Can only happen if the root is unsupported.
    cur_node_ptr->is_unsupported = true;
}

}  // namespace algos::hymd::lattice

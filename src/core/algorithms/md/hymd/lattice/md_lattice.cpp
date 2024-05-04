#include "algorithms/md/hymd/lattice/md_lattice.h"

#include <algorithm>
#include <cassert>

#include "algorithms/md/hymd/lattice/spec_generalization_checker.h"
#include "algorithms/md/hymd/lattice/total_generalization_checker.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/utility/set_for_scope.h"
#include "util/erase_if_replace.h"

namespace {
using namespace algos::hymd::lattice;
using namespace algos::hymd;
using model::md::DecisionBoundary, model::Index;
using MdGenChecker = TotalGeneralizationChecker<MdNode>;
using MdSpecGenChecker = SpecGeneralizationChecker<MdNode>;

bool NotEmpty(DecisionBoundaryVector const& rhs_bounds) {
    auto not_lowest = [](DecisionBoundary bound) { return bound != kLowestBound; };
    return std::any_of(rhs_bounds.begin(), rhs_bounds.end(), not_lowest);
}
}  // namespace

namespace algos::hymd::lattice {

// TODO: remove recursion
MdLattice::MdLattice(std::size_t column_matches_size, SingleLevelFunc single_level_func,
                     std::vector<std::vector<model::md::DecisionBoundary>> const& lhs_bounds,
                     bool prune_nondisjoint)
    : column_matches_size_(column_matches_size),
      md_root_(DecisionBoundaryVector(column_matches_size_, 1.0)),
      support_root_(column_matches_size_),
      get_single_level_(std::move(single_level_func)),
      lhs_bounds_(&lhs_bounds),
      prune_nondisjoint_(prune_nondisjoint) {}

std::optional<DecisionBoundary> MdLattice::SpecializeOneLhs(Index col_match_index,
                                                            DecisionBoundary lhs_bound) const {
    std::vector<DecisionBoundary> const& decision_bounds = (*lhs_bounds_)[col_match_index];
    auto end_bounds = decision_bounds.end();
    auto upper = std::upper_bound(decision_bounds.begin(), end_bounds, lhs_bound);
    if (upper == end_bounds) return std::nullopt;
    return *upper;
}

void MdLattice::Specialize(MdLhs const& lhs, SimilarityVector const& specialize_past,
                           Rhss const& rhss) {
    auto specialize_all_lhs = [this, &lhs, it_begin = rhss.begin(), it_end = rhss.end(),
                               &specialize_past](auto handle_same_lhs_as_rhs) {
        for (Index lhs_spec_index = 0; lhs_spec_index < column_matches_size_; ++lhs_spec_index) {
            std::optional<DecisionBoundary> const specialized_lhs_bound =
                    SpecializeOneLhs(lhs_spec_index, specialize_past[lhs_spec_index]);
            if (!specialized_lhs_bound.has_value()) continue;

            LhsSpecialization lhs_specialization{lhs, {lhs_spec_index, *specialized_lhs_bound}};
            if (IsUnsupported(lhs_specialization)) continue;

            for (auto it = it_begin; it != it_end; ++it) {
                if (it->index == lhs_spec_index) {
                    handle_same_lhs_as_rhs(lhs_specialization, *it);
                    for (++it; it != it_end; ++it) {
                        AddIfMinimal({lhs_specialization, *it});
                    }
                    break;
                }
                AddIfMinimal({lhs_specialization, *it});
            }
        }
    };
    if (prune_nondisjoint_) {
        specialize_all_lhs([](...) {});
    } else {
        specialize_all_lhs([this](LhsSpecialization const& lhs_spec, MdElement rhs) {
            if (rhs.decision_boundary > lhs_spec.specialized.decision_boundary) {
                AddIfMinimal({lhs_spec, rhs});
            }
        });
    }
}

void MdLattice::Specialize(MdLhs const& lhs, Rhss const& rhss) {
    Specialize(lhs, lhs.GetValues(), rhss);
}

void MdLattice::MdRefiner::Refine() {
    for (auto upd_iter = invalidated_.UpdateIterBegin(), upd_end = invalidated_.UpdateIterEnd();
         upd_iter != upd_end; ++upd_iter) {
        auto const& [rhs_index, new_bound] = *upd_iter;
        DecisionBoundary& md_rhs_bound_ref = (*node_info_.rhs_bounds)[rhs_index];
        md_rhs_bound_ref = kLowestBound;
        // trivial
        if (new_bound == GetLhs()[rhs_index]) continue;
        // not minimal
        if (lattice_->HasGeneralization({GetLhs(), *upd_iter})) continue;
        md_rhs_bound_ref = new_bound;
    }
    lattice_->Specialize(GetLhs(), *sim_, invalidated_.GetInvalidated());
}

void MdLattice::TryAddRefiner(std::vector<MdRefiner>& found, DecisionBoundaryVector& rhs,
                              SimilarityVector const& similarity_vector,
                              MdLhs const& cur_node_lhs) {
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
        found.emplace_back(this, &similarity_vector, MdLatticeNodeInfo{cur_node_lhs, &rhs},
                           std::move(invalidated));
        break;
    }
}

void MdLattice::CollectRefinersForViolated(MdNode& cur_node, std::vector<MdRefiner>& found,
                                           MdLhs& cur_node_lhs,
                                           SimilarityVector const& similarity_vector,
                                           Index cur_node_index) {
    TryAddRefiner(found, cur_node.rhs_bounds, similarity_vector, cur_node_lhs);

    auto collect = [&](MdBoundMap& bound_map, model::Index child_array_index) {
        Index const next_node_index = cur_node_index + child_array_index;
        DecisionBoundary& cur_lhs_bound = cur_node_lhs[next_node_index];
        DecisionBoundary const sim_vec_sim = similarity_vector[next_node_index];
        for (auto& [generalization_bound, node] : bound_map) {
            if (generalization_bound > sim_vec_sim) break;
            cur_lhs_bound = generalization_bound;
            CollectRefinersForViolated(node, found, cur_node_lhs, similarity_vector,
                                       next_node_index + 1);
        }
        cur_lhs_bound = kLowestBound;
    };
    cur_node.ForEachNonEmpty(collect);
}

auto MdLattice::CollectRefinersForViolated(SimilarityVector const& similarity_vector)
        -> std::vector<MdRefiner> {
    std::vector<MdRefiner> found;
    MdLhs current_lhs(column_matches_size_);
    CollectRefinersForViolated(md_root_, found, current_lhs, similarity_vector, 0);
    // TODO: traverse support trie simultaneously.
    util::EraseIfReplace(found, [this](MdRefiner& refiner) {
        bool const unsupported = IsUnsupported(refiner.GetLhs());
        if (unsupported) refiner.ZeroRhs();
        return unsupported;
    });
    return found;
}

bool MdLattice::IsUnsupported(MdLhs const& lhs) const {
    return TotalGeneralizationChecker<SupportNode>{lhs}.HasGeneralization(support_root_);
}

void MdLattice::MdVerificationMessenger::MarkUnsupported() {
    // TODO: specializations can be removed from the MD lattice. If not worth it, removing just
    // this node and its children should be cheap. Though, destructors also take time.

    // This matters. Violation search can find a node with a specialized LHS but higher RHS bound,
    // leading to extra work (though no influence on correctness, as MDs with unsupported LHSs are
    // filtered out).
    ZeroRhs();

    lattice_->MarkUnsupported(GetLhs());
}

void MdLattice::MdVerificationMessenger::LowerAndSpecialize(
        utility::InvalidatedRhss const& invalidated) {
    for (auto upd_iter = invalidated.UpdateIterBegin(), upd_end = invalidated.UpdateIterEnd();
         upd_iter != upd_end; ++upd_iter) {
        auto const& [rhs_index, new_bound] = *upd_iter;
        GetRhs()[rhs_index] = new_bound;
    }
    lattice_->Specialize(GetLhs(), invalidated.GetInvalidated());
}

void MdLattice::AddNewMinimal(MdNode& cur_node, MdSpecialization const& md, Index cur_node_index) {
    assert(!NotEmpty(cur_node.rhs_bounds));
    assert(cur_node_index > md.lhs_specialization.specialized.index);
    auto const& [rhs_index, rhs_bound] = md.rhs;
    auto set_bound = [&](MdNode* node) { node->rhs_bounds[rhs_index] = rhs_bound; };
    AddUnchecked(&cur_node, md.lhs_specialization.old_lhs, cur_node_index, set_bound);
    UpdateMaxLevel(md.lhs_specialization);
}

void MdLattice::UpdateMaxLevel(LhsSpecialization const& lhs) {
    std::size_t level = 0;
    auto const& [spec_index, spec_bound] = lhs.specialized;
    MdLhs const& old_lhs = lhs.old_lhs;
    for (Index i = 0; i != spec_index; ++i) {
        DecisionBoundary const cur_bound = old_lhs[i];
        if (cur_bound == kLowestBound) continue;
        level += get_single_level_(cur_bound, i);
    }
    assert(spec_bound != kLowestBound);
    level += get_single_level_(spec_bound, spec_index);
    for (Index i = spec_index + 1; i != column_matches_size_; ++i) {
        DecisionBoundary const cur_bound = old_lhs[i];
        if (cur_bound == kLowestBound) continue;
        level += get_single_level_(cur_bound, i);
    }
    if (level > max_level_) max_level_ = level;
}

class MdLattice::GeneralizationHelper {
private:
    MdElement const rhs_;
    MdNode* node_;
    DecisionBoundary* cur_node_rhs_ptr_;
    MdGenChecker const& gen_checker_;

public:
    GeneralizationHelper(MdElement rhs, MdNode& root, MdGenChecker const& gen_checker) noexcept
        : rhs_(rhs),
          node_(&root),
          cur_node_rhs_ptr_(&node_->rhs_bounds[rhs_.index]),
          gen_checker_(gen_checker) {}

    bool SetAndCheck(MdNode* node_ptr) noexcept {
        if (!node_ptr) return true;
        node_ = node_ptr;
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

    MdGenChecker const& GetTotalChecker() const noexcept {
        return gen_checker_;
    }
};

// Note: writing this in AddIfMinimal with gotos seems to be faster.
auto MdLattice::TryGetNextNode(MdSpecialization const& md, GeneralizationHelper& helper,
                               Index cur_node_index, Index const next_node_index,
                               DecisionBoundary const next_lhs_bound) -> MdNode* {
    Index const fol_index = next_node_index + 1;
    Index const child_array_index = next_node_index - cur_node_index;
    MdNode& cur_node = helper.CurNode();
    auto [boundary_mapping, is_first_arr] = cur_node.TryEmplaceChild(child_array_index);
    std::size_t const next_child_array_size = cur_node.GetChildArraySize(child_array_index);
    if (is_first_arr) [[unlikely]] {
        MdNode& new_node =
                boundary_mapping
                        .try_emplace(next_lhs_bound, column_matches_size_, next_child_array_size)
                        .first->second;
        AddNewMinimal(new_node, md, fol_index);
        return nullptr;
    }
    auto it = boundary_mapping.begin();
    MdGenChecker total_checker = helper.GetTotalChecker();
    for (auto end_it = boundary_mapping.end(); it != end_it; ++it) {
        auto const& [generalization_bound, node] = *it;
        if (generalization_bound > next_lhs_bound) break;
        if (generalization_bound == next_lhs_bound) return &it->second;
        if (total_checker.HasGeneralization(node, fol_index)) return nullptr;
    }
    using std::forward_as_tuple;
    MdNode& new_node =
            boundary_mapping
                    .emplace_hint(it, std::piecewise_construct, forward_as_tuple(next_lhs_bound),
                                  forward_as_tuple(column_matches_size_, next_child_array_size))
                    ->second;
    AddNewMinimal(new_node, md, fol_index);
    return nullptr;
}

void MdLattice::AddIfMinimal(MdSpecialization const& md) {
    MdSpecGenChecker gen_checker{md};
    MdGenChecker const& total_checker = gen_checker.GetTotalChecker();
    auto helper = GeneralizationHelper(md.rhs, md_root_, total_checker);
    Index cur_node_index = 0;
    auto const& [spec_index, spec_bound] = md.lhs_specialization.specialized;
    MdLhs const& old_lhs = md.lhs_specialization.old_lhs;
    auto try_set_next = [&](auto... args) {
        return helper.SetAndCheck(TryGetNextNode(md, helper, cur_node_index, args...));
    };
    for (MdElement result = old_lhs.FindNextNonZero(cur_node_index); result.index < spec_index;
         result = old_lhs.FindNextNonZero(cur_node_index)) {
        auto const& [next_node_index, next_lhs_bound] = result;
        Index const child_array_index = next_node_index - cur_node_index;
        Index const fol_index = result.index + 1;
        if (gen_checker.HasGeneralizationInChildren(helper.CurNode(), cur_node_index, fol_index))
            return;
        cur_node_index = fol_index;
        assert(helper.Children()[child_array_index].has_value());
        MdBoundMap& bound_map = *helper.Children()[child_array_index];
        assert(bound_map.find(next_lhs_bound) != bound_map.end());
        auto it = bound_map.begin();
        for (; it->first != next_lhs_bound; ++it) {
            if (gen_checker.HasGeneralization(it->second, cur_node_index)) return;
        }
        helper.SetAndCheck(&it->second);
    }
    if (try_set_next(spec_index, spec_bound)) return;

    cur_node_index = spec_index + 1;
    for (MdElement element = old_lhs.FindNextNonZero(cur_node_index); old_lhs.IsNotEnd(element);
         element = old_lhs.FindNextNonZero(cur_node_index)) {
        auto const& [next_node_index, next_lhs_bound] = element;
        Index const fol_index = next_node_index + 1;
        if (total_checker.HasGeneralizationInChildren(helper.CurNode(), cur_node_index, fol_index))
            return;
        if (try_set_next(next_node_index, next_lhs_bound)) return;
        cur_node_index = fol_index;
    }
    // Note: Metanome implemented this incorrectly, potentially missing out on recommendations.
    helper.SetBoundOnCurrent();
}

void MdLattice::RaiseInterestingnessBounds(
        MdNode const& cur_node, MdLhs const& lhs,
        std::vector<DecisionBoundary>& cur_interestingness_bounds, Index const cur_node_index,
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

    for (MdElement element = lhs.FindNextNonZero(cur_node_index); lhs.IsNotEnd(element);
         element = lhs.FindNextNonZero(element.index + 1)) {
        auto const& [next_node_index, generalization_bound_limit] = element;
        Index const child_array_index = next_node_index - cur_node_index;
        MdOptionalChild const& optional_child = cur_node.children[child_array_index];
        if (!optional_child.has_value()) continue;
        for (auto const& [generalization_bound, node] : *optional_child) {
            if (generalization_bound > generalization_bound_limit) break;
            RaiseInterestingnessBounds(node, lhs, cur_interestingness_bounds, next_node_index + 1,
                                       indices);
        }
    }
}

std::vector<DecisionBoundary> MdLattice::GetRhsInterestingnessBounds(
        MdLhs const& lhs, std::vector<Index> const& indices) const {
    std::vector<DecisionBoundary> interestingness_bounds;
    interestingness_bounds.reserve(indices.size());
    for (Index index : indices) {
        DecisionBoundary const lhs_bound = lhs[index];
        assert(lhs_bound != 1.0);
        interestingness_bounds.push_back(lhs_bound);
    }
    RaiseInterestingnessBounds(md_root_, lhs, interestingness_bounds, 0, indices);
    return interestingness_bounds;
}

bool MdLattice::HasGeneralization(Md const& md) const {
    return MdGenChecker{md}.HasGeneralization(md_root_);
}

void MdLattice::GetLevel(MdNode& cur_node, std::vector<MdVerificationMessenger>& collected,
                         MdLhs& cur_node_lhs, Index const cur_node_index,
                         std::size_t const level_left) {
    DecisionBoundaryVector& rhs_bounds = cur_node.rhs_bounds;
    if (level_left == 0) {
        if (NotEmpty(rhs_bounds))
            collected.emplace_back(this, MdLatticeNodeInfo{cur_node_lhs, &rhs_bounds});
        return;
    }
    auto collect = [&](MdBoundMap& bound_map, model::Index child_array_index) {
        Index const next_node_index = cur_node_index + child_array_index;
        DecisionBoundary& next_lhs_bound = cur_node_lhs[next_node_index];
        for (auto& [boundary, node] : bound_map) {
            std::size_t const single = get_single_level_(next_node_index, boundary);
            if (single > level_left) break;
            next_lhs_bound = boundary;
            GetLevel(node, collected, cur_node_lhs, next_node_index + 1, level_left - single);
        }
        next_lhs_bound = kLowestBound;
    };
    cur_node.ForEachNonEmpty(collect);
}

auto MdLattice::GetLevel(std::size_t const level) -> std::vector<MdVerificationMessenger> {
    std::vector<MdVerificationMessenger> collected;
    MdLhs current_lhs(column_matches_size_);
    GetLevel(md_root_, collected, current_lhs, 0, level);
    // TODO: traverse support trie simultaneously.
    util::EraseIfReplace(collected, [this](MdVerificationMessenger& messenger) {
        bool is_unsupported = IsUnsupported(messenger.GetLhs());
        if (is_unsupported) messenger.ZeroRhs();
        return is_unsupported;
    });
    return collected;
}

void MdLattice::GetAll(MdNode& cur_node, std::vector<MdLatticeNodeInfo>& collected,
                       MdLhs& cur_node_lhs, Index const this_node_index) {
    DecisionBoundaryVector& rhs_bounds = cur_node.rhs_bounds;
    if (NotEmpty(rhs_bounds)) collected.emplace_back(cur_node_lhs, &rhs_bounds);
    auto collect = [&](MdBoundMap& bound_map, model::Index child_array_index) {
        Index const next_node_index = this_node_index + child_array_index;
        DecisionBoundary& next_lhs_bound = cur_node_lhs[next_node_index];
        for (auto& [boundary, node] : bound_map) {
            next_lhs_bound = boundary;
            GetAll(node, collected, cur_node_lhs, next_node_index + 1);
        }
        next_lhs_bound = kLowestBound;
    };
    cur_node.ForEachNonEmpty(collect);
}

std::vector<MdLatticeNodeInfo> MdLattice::GetAll() {
    std::vector<MdLatticeNodeInfo> collected;
    MdLhs current_lhs(column_matches_size_);
    GetAll(md_root_, collected, current_lhs, 0);
    assert(std::none_of(
            collected.begin(), collected.end(),
            [this](MdLatticeNodeInfo const& node_info) { return IsUnsupported(node_info.lhs); }));
    return collected;
}

bool MdLattice::IsUnsupported(LhsSpecialization const& lhs_spec) const {
    return SpecGeneralizationChecker<SupportNode>{lhs_spec}.HasGeneralization(support_root_, 0);
}

void MdLattice::MarkNewLhs(SupportNode& cur_node, MdLhs const& lhs, Index cur_node_index) {
    AddUnchecked(&cur_node, lhs, cur_node_index, SetUnsupAction());
}

void MdLattice::MarkUnsupported(MdLhs const& lhs) {
    auto mark_new = [this](auto&&... args) { MarkNewLhs(std::forward<decltype(args)>(args)...); };
    CheckedAdd(&support_root_, lhs, lhs, mark_new, SetUnsupAction());
}

}  // namespace algos::hymd::lattice

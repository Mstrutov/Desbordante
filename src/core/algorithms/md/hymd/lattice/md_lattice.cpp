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

void MdLattice::Specialize(MdLhs const& lhs, Rhss const& rhss, auto get_higher_lhs_bound,
                           auto get_higher_other_bound) {
    Index lhs_spec_index = 0;
    auto rhs_begin = rhss.begin(), rhs_end = rhss.end();
    auto lhs_iter = lhs.begin(), lhs_end = lhs.end();
    auto add_all_rhs = [&](LhsSpecialization const& lhs_spec) {
        for (auto rhs_it = rhs_begin; rhs_it != rhs_end; ++rhs_it) {
            if (rhs_it->index == lhs_spec_index) {
                if (!prune_nondisjoint_) {
                    DecisionBoundary const specialized_lhs_bound =
                            lhs_spec.specialization_data.new_child.decision_boundary;
                    if (specialized_lhs_bound < rhs_it->decision_boundary)
                        AddIfMinimal({lhs_spec, *rhs_it});
                }
                for (++rhs_it; rhs_it != rhs_end; ++rhs_it) {
                    AddIfMinimal({lhs_spec, *rhs_it});
                }
                return;
            }
            AddIfMinimal({lhs_spec, *rhs_it});
        }
    };
    auto specialize_element = [&](Index spec_child_index, DecisionBoundary spec_past) {
        std::optional<DecisionBoundary> const specialized_lhs_bound =
                SpecializeOneLhs(lhs_spec_index, spec_past);
        if (!specialized_lhs_bound.has_value()) return;
        LhsSpecialization lhs_spec{lhs, {lhs_iter, {spec_child_index, *specialized_lhs_bound}}};
        if (IsUnsupported(lhs_spec)) return;
        add_all_rhs(lhs_spec);
    };
    for (; lhs_iter != lhs_end; ++lhs_iter, ++lhs_spec_index) {
        auto const& [child_array_index, bound] = *lhs_iter;
        for (Index spec_child_index = 0; spec_child_index != child_array_index;
             ++spec_child_index, ++lhs_spec_index) {
            specialize_element(spec_child_index, get_higher_other_bound(lhs_spec_index));
        }
        specialize_element(child_array_index, get_higher_lhs_bound(lhs_spec_index, bound));
    }
    for (Index spec_child_index = 0; lhs_spec_index != column_matches_size_;
         ++lhs_spec_index, ++spec_child_index) {
        specialize_element(spec_child_index, get_higher_other_bound(lhs_spec_index));
    };
}

void MdLattice::Specialize(MdLhs const& lhs, SimilarityVector const& specialize_past,
                           Rhss const& rhss) {
    auto get_sim_vec_bound = [&](Index index, ...) { return specialize_past[index]; };
    Specialize(lhs, rhss, get_sim_vec_bound, get_sim_vec_bound);
}

void MdLattice::Specialize(MdLhs const& lhs, Rhss const& rhss) {
    auto get_lowest = [&](...) { return kLowestBound; };
    auto get_lhs_bound = [](Index, DecisionBoundary bound) { return bound; };
    Specialize(lhs, rhss, get_lhs_bound, get_lowest);
}

void MdLattice::MdRefiner::Refine() {
    for (auto upd_iter = invalidated_.UpdateIterBegin(), upd_end = invalidated_.UpdateIterEnd();
         upd_iter != upd_end; ++upd_iter) {
        auto const& [rhs_index, new_bound] = *upd_iter;
        DecisionBoundary& md_rhs_bound_ref = (*node_info_.rhs_bounds)[rhs_index];
        md_rhs_bound_ref = kLowestBound;
        // trivial
        if (new_bound == kLowestBound) continue;
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
    Index rhs_index = 0;
    Index cur_lhs_index = 0;
    auto try_push_no_match_classifier = [&]() {
        DecisionBoundary sim_bound = similarity_vector[rhs_index];
        DecisionBoundary rhs_bound = rhs[rhs_index];
        if (sim_bound < rhs_bound) {
            invalidated.PushBack({rhs_index, rhs_bound}, sim_bound);
        }
    };
    for (auto const& [child_index, lhs_bound] : cur_node_lhs) {
        cur_lhs_index += child_index;
        for (; rhs_index != cur_lhs_index; ++rhs_index) {
            try_push_no_match_classifier();
        }
        assert(rhs_index < column_matches_size_);
        DecisionBoundary const sim_bound = similarity_vector[rhs_index];
        DecisionBoundary const rhs_bound = rhs[rhs_index];
        if (sim_bound < rhs_bound) {
            MdElement invalid{rhs_index, rhs_bound};
            if (sim_bound == lhs_bound) {
                invalidated.PushBack(invalid, kLowestBound);
            } else {
                assert(sim_bound > lhs_bound);
                invalidated.PushBack(invalid, sim_bound);
            }
        }
        ++cur_lhs_index;
    }
    for (; rhs_index != column_matches_size_; ++rhs_index) {
        try_push_no_match_classifier();
    }
    if (invalidated.IsEmpty()) return;
    found.emplace_back(this, &similarity_vector, MdLatticeNodeInfo{cur_node_lhs, &rhs},
                       std::move(invalidated));
}

void MdLattice::CollectRefinersForViolated(MdNode& cur_node, std::vector<MdRefiner>& found,
                                           MdLhs& cur_node_lhs,
                                           SimilarityVector const& similarity_vector,
                                           Index cur_node_index) {
    TryAddRefiner(found, cur_node.rhs_bounds, similarity_vector, cur_node_lhs);

    auto collect = [&](MdBoundMap& bound_map, model::Index child_array_index) {
        Index const next_node_index = cur_node_index + child_array_index;
        DecisionBoundary& cur_lhs_bound = cur_node_lhs.AddNext(child_array_index);
        DecisionBoundary const sim_vec_sim = similarity_vector[next_node_index];
        for (auto& [generalization_bound, node] : bound_map) {
            if (generalization_bound > sim_vec_sim) break;
            cur_lhs_bound = generalization_bound;
            CollectRefinersForViolated(node, found, cur_node_lhs, similarity_vector,
                                       next_node_index + 1);
        }
        cur_node_lhs.RemoveLast();
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
    // assert(cur_node_index > md.lhs_specialization.specialized.index);
    auto const& [rhs_index, rhs_bound] = md.rhs;
    auto set_bound = [&](MdNode* node) { node->rhs_bounds[rhs_index] = rhs_bound; };
    AddUnchecked(&cur_node, md.lhs_specialization.old_lhs, cur_node_index, set_bound);
    UpdateMaxLevel(md.lhs_specialization);
}

void MdLattice::AddNewMinimal(MdNode& cur_node, MdSpecialization const& md,
                              MdLhs::iterator cur_node_iter) {
    assert(!NotEmpty(cur_node.rhs_bounds));
    assert(cur_node_iter >= md.lhs_specialization.specialization_data.spec_before);
    auto const& [rhs_index, rhs_bound] = md.rhs;
    auto set_bound = [&](MdNode* node) { node->rhs_bounds[rhs_index] = rhs_bound; };
    AddUnchecked(&cur_node, md.lhs_specialization.old_lhs, cur_node_iter, set_bound);
    UpdateMaxLevel(md.lhs_specialization);
}

void MdLattice::UpdateMaxLevel(LhsSpecialization const& lhs) {
    std::size_t level = 0;
    auto const& [spec_child_array_index, spec_bound] = lhs.specialization_data.new_child;
    MdLhs const& old_lhs = lhs.old_lhs;
    MdLhs::iterator hint_iter = lhs.specialization_data.spec_before;
    Index cur_col_match_index = 0;
    MdLhs::iterator lhs_iter = old_lhs.begin();
    auto add_level = [&]() {
        auto const& [index_delta, bound] = *lhs_iter;
        cur_col_match_index += index_delta;
        level += get_single_level_(bound, cur_col_match_index);
        ++cur_col_match_index;
    };
    for (; lhs_iter != hint_iter; ++lhs_iter) {
        add_level();
    }
    level += get_single_level_(spec_bound, cur_col_match_index + spec_child_array_index);
    MdLhs::iterator lhs_end = old_lhs.end();
    if (lhs_iter != lhs_end) {
        if (spec_child_array_index != lhs_iter->child_array_index) add_level();
        for (++lhs_iter; lhs_iter != lhs_end; ++lhs_iter) {
            add_level();
        }
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
MdNode* MdLattice::TryGetNextNode(MdSpecialization const&, GeneralizationHelper& helper,
                                  Index const child_array_index, auto new_minimal_action,
                                  DecisionBoundary const next_lhs_bound, MdLhs::iterator iter,
                                  std::size_t gen_check_offset) {
    MdNode& cur_node = helper.CurNode();
    auto [boundary_mapping, is_first_arr] = cur_node.TryEmplaceChild(child_array_index);
    std::size_t const next_child_array_size = cur_node.GetChildArraySize(child_array_index);
    if (is_first_arr) [[unlikely]] {
        MdNode& new_node =
                boundary_mapping
                        .try_emplace(next_lhs_bound, column_matches_size_, next_child_array_size)
                        .first->second;
        new_minimal_action(new_node);
        return nullptr;
    }
    auto it = boundary_mapping.begin();
    MdGenChecker total_checker = helper.GetTotalChecker();
    for (auto end_it = boundary_mapping.end(); it != end_it; ++it) {
        auto const& [generalization_bound, node] = *it;
        if (generalization_bound > next_lhs_bound) break;
        if (generalization_bound == next_lhs_bound) return &it->second;
        if (total_checker.HasGeneralization(node, iter, gen_check_offset)) return nullptr;
    }
    using std::forward_as_tuple;
    MdNode& new_node =
            boundary_mapping
                    .emplace_hint(it, std::piecewise_construct, forward_as_tuple(next_lhs_bound),
                                  forward_as_tuple(column_matches_size_, next_child_array_size))
                    ->second;
    new_minimal_action(new_node);
    return nullptr;
}

void MdLattice::AddIfMinimal(MdSpecialization const& md) {
    MdSpecGenChecker gen_checker{md};
    MdGenChecker const& total_checker = gen_checker.GetTotalChecker();
    auto helper = GeneralizationHelper(md.rhs, md_root_, total_checker);
    auto const& [spec_child_array_index, spec_bound] =
            md.lhs_specialization.specialization_data.new_child;
    MdLhs::iterator spec_iter = md.lhs_specialization.specialization_data.spec_before;
    MdLhs const& old_lhs = md.lhs_specialization.old_lhs;
    Index const spec_index = old_lhs.GetColumnMatchIndex(spec_iter, spec_child_array_index);
    MdLhs::iterator next_lhs_iter = old_lhs.begin();
    while (next_lhs_iter != spec_iter) {
        auto const& [child_array_index, next_lhs_bound] = *next_lhs_iter;
        ++next_lhs_iter;
        if (gen_checker.HasGeneralizationInChildren(helper.CurNode(), next_lhs_iter,
                                                    child_array_index + 1))
            return;
        assert(helper.Children()[child_array_index].has_value());
        MdBoundMap& bound_map = *helper.Children()[child_array_index];
        assert(bound_map.find(next_lhs_bound) != bound_map.end());
        auto it = bound_map.begin();
        for (; it->first != next_lhs_bound; ++it) {
            if (gen_checker.HasGeneralization(it->second, next_lhs_iter)) return;
        }
        helper.SetAndCheck(&it->second);
    }
    auto try_set_next = [&](auto... args) {
        return helper.SetAndCheck(TryGetNextNode(md, helper, args...));
    };

    MdLhs::iterator lhs_end = old_lhs.end();
    if (next_lhs_iter == lhs_end) {  // Append
        auto new_minimal_action = [&](MdNode& node) { AddNewMinimal(node, md, spec_index + 1); };
        if (try_set_next(spec_child_array_index, new_minimal_action, spec_bound, next_lhs_iter))
            return;
    } else if (next_lhs_iter->child_array_index == spec_child_array_index) {  // Replace
        ++next_lhs_iter;
        auto new_minimal_action = [&](MdNode& node) { AddNewMinimal(node, md, spec_index + 1); };
        if (try_set_next(spec_child_array_index, new_minimal_action, spec_bound, next_lhs_iter))
            return;
    } else {  // Insert
        auto const& [old_child_array_index, next_lhs_bound] = *next_lhs_iter;
        std::size_t const offset = -(spec_child_array_index + 1);
        auto new_minimal_action = [&](MdNode& node) { AddNewMinimal(node, md, spec_index + 1); };
        if (try_set_next(spec_child_array_index, new_minimal_action, spec_bound, next_lhs_iter,
                         offset))
            return;
        if (total_checker.HasGeneralizationInChildren(helper.CurNode(), next_lhs_iter, offset))
            return;
        std::size_t const fol_spec_child_index = old_child_array_index + offset;
        model::Index next_node_index = old_lhs.ToIndex(next_lhs_iter);
        ++next_lhs_iter;
        auto new_minimal_action2 = [&](MdNode& node) {
            AddNewMinimal(node, md, next_node_index + 1);
        };
        if (try_set_next(fol_spec_child_index, new_minimal_action2, next_lhs_bound, next_lhs_iter))
            return;
    }
    while (next_lhs_iter != lhs_end) {
        auto const& [child_array_index, next_lhs_bound] = *next_lhs_iter;
        model::Index next_node_index = old_lhs.ToIndex(next_lhs_iter);
        auto new_minimal_action2 = [&](MdNode& node) {
            AddNewMinimal(node, md, next_node_index + 1);
        };
        ++next_lhs_iter;
        if (total_checker.HasGeneralizationInChildren(helper.CurNode(), next_lhs_iter,
                                                      child_array_index + 1))
            return;
        if (try_set_next(child_array_index, new_minimal_action2, next_lhs_bound, next_lhs_iter))
            return;
    }
    // Note: Metanome implemented this incorrectly, potentially missing out on recommendations.
    helper.SetBoundOnCurrent();
}

void MdLattice::RaiseInterestingnessBounds(
        MdNode const& cur_node, MdLhs const& lhs,
        std::vector<DecisionBoundary>& cur_interestingness_bounds, MdLhs::iterator cur_lhs_iter,
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

    Index child_array_index = 0;
    for (MdLhs::iterator end = lhs.end(); cur_lhs_iter != end;
         ++cur_lhs_iter, ++child_array_index) {
        auto const& [offset, generalization_bound_limit] = *cur_lhs_iter;
        child_array_index += offset;
        MdOptionalChild const& optional_child = cur_node.children[child_array_index];
        if (!optional_child.has_value()) continue;
        for (auto const& [generalization_bound, node] : *optional_child) {
            if (generalization_bound > generalization_bound_limit) break;
            RaiseInterestingnessBounds(node, lhs, cur_interestingness_bounds, cur_lhs_iter + 1,
                                       indices);
        }
    }
}

std::vector<DecisionBoundary> MdLattice::GetRhsInterestingnessBounds(
        MdLhs const& lhs, std::vector<Index> const& indices) const {
    std::vector<DecisionBoundary> interestingness_bounds;
    std::size_t indices_size = indices.size();
    if (prune_nondisjoint_) {
        interestingness_bounds.assign(indices_size, kLowestBound);
    } else {
        interestingness_bounds.reserve(indices_size);
        assert(std::is_sorted(indices.begin(), indices.end()));
        auto fill_bounds = [&]() {
            auto index_it = indices.begin(), index_end = indices.end();
            Index cur_index = 0;
            assert(!indices.empty());
            for (auto const& [child_index, bound] : lhs) {
                cur_index += child_index;
                Index index;
                while ((index = *index_it) < cur_index) {
                    interestingness_bounds.push_back(kLowestBound);
                    if (++index_it == index_end) return;
                }
                if (cur_index == index) {
                    interestingness_bounds.push_back(bound);
                }
                ++cur_index;
            }
            for (; index_it != index_end; ++index_it) {
                interestingness_bounds.push_back(kLowestBound);
            }
        };
        fill_bounds();
    }
    RaiseInterestingnessBounds(md_root_, lhs, interestingness_bounds, lhs.begin(), indices);
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
        DecisionBoundary& next_lhs_bound = cur_node_lhs.AddNext(child_array_index);
        for (auto& [boundary, node] : bound_map) {
            std::size_t const single = get_single_level_(next_node_index, boundary);
            if (single > level_left) break;
            next_lhs_bound = boundary;
            GetLevel(node, collected, cur_node_lhs, next_node_index + 1, level_left - single);
        }
        cur_node_lhs.RemoveLast();
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
        DecisionBoundary& next_lhs_bound = cur_node_lhs.AddNext(child_array_index);
        for (auto& [boundary, node] : bound_map) {
            next_lhs_bound = boundary;
            GetAll(node, collected, cur_node_lhs, next_node_index + 1);
        }
        cur_node_lhs.RemoveLast();
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
    return SpecGeneralizationChecker<SupportNode>{lhs_spec}.HasGeneralization(support_root_);
}

void MdLattice::MarkNewLhs(SupportNode& cur_node, MdLhs const& lhs, MdLhs::iterator cur_lhs_iter) {
    AddUnchecked(&cur_node, lhs, cur_lhs_iter, SetUnsupAction());
}

void MdLattice::MarkUnsupported(MdLhs const& lhs) {
    auto mark_new = [this](auto&&... args) { MarkNewLhs(std::forward<decltype(args)>(args)...); };
    CheckedAdd(&support_root_, lhs, lhs, mark_new, SetUnsupAction());
}

}  // namespace algos::hymd::lattice

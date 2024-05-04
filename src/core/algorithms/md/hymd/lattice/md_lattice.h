#pragma once

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lattice/md.h"
#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/md_node.h"
#include "algorithms/md/hymd/lattice/node_base.h"
#include "algorithms/md/hymd/lattice/single_level_func.h"
#include "algorithms/md/hymd/lattice/support_node.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/md_lhs.h"
#include "algorithms/md/hymd/rhss.h"
#include "algorithms/md/hymd/similarity_vector.h"
#include "algorithms/md/hymd/utility/invalidated_rhss.h"
#include "model/index.h"

namespace algos::hymd::lattice {

class MdLattice {
private:
    class GeneralizationHelper;

    using MdBoundMap = MdNode::BoundMap;
    using MdOptionalChild = MdNode::OptionalChild;
    using MdNodeChildren = MdNode::Children;

    using SupportBoundMap = SupportNode::BoundMap;
    using SupportOptionalChild = SupportNode::OptionalChild;
    using SupportNodeChildren = SupportNode::Children;

public:
    class MdRefiner {
        MdLattice* lattice_;
        SimilarityVector const* sim_;
        MdLatticeNodeInfo node_info_;
        utility::InvalidatedRhss invalidated_;

    public:
        MdRefiner(MdLattice* lattice, SimilarityVector const* sim, MdLatticeNodeInfo node_info,
                  utility::InvalidatedRhss invalidated)
            : lattice_(lattice),
              sim_(sim),
              node_info_(std::move(node_info)),
              invalidated_(std::move(invalidated)) {}

        MdLhs const& GetLhs() const {
            return node_info_.lhs;
        }

        void ZeroRhs() {
            node_info_.ZeroRhs();
        }

        void Refine();
    };

    class MdVerificationMessenger {
        MdLattice* lattice_;
        MdLatticeNodeInfo node_info_;

    public:
        MdVerificationMessenger(MdLattice* lattice, MdLatticeNodeInfo node_info)
            : lattice_(lattice), node_info_(std::move(node_info)) {}

        MdLhs const& GetLhs() const {
            return node_info_.lhs;
        }

        DecisionBoundaryVector& GetRhs() {
            return *node_info_.rhs_bounds;
        }

        void MarkUnsupported();

        void ZeroRhs() {
            node_info_.ZeroRhs();
        }

        void LowerAndSpecialize(utility::InvalidatedRhss const& invalidated);
    };

private:
    std::size_t max_level_ = 0;
    std::size_t const column_matches_size_;
    MdNode md_root_;
    SupportNode support_root_;
    // Is there a way to define a level in such a way that one cannot use each decision boundary
    // independently to determine an MD's level but the lattice traversal algorithms still works?
    SingleLevelFunc const get_single_level_;
    std::vector<std::vector<model::md::DecisionBoundary>> const* const lhs_bounds_;
    bool const prune_nondisjoint_;

    [[nodiscard]] bool HasGeneralization(Md const& md) const;

    void GetLevel(MdNode& cur_node, std::vector<MdVerificationMessenger>& collected,
                  MdLhs& cur_node_lhs, model::Index cur_node_index, std::size_t level_left);

    void RaiseInterestingnessBounds(
            MdNode const& cur_node, MdLhs const& lhs,
            std::vector<model::md::DecisionBoundary>& cur_interestingness_bounds,
            model::Index cur_node_index, std::vector<model::Index> const& indices) const;

    void TryAddRefiner(std::vector<MdRefiner>& found, DecisionBoundaryVector& rhs,
                       SimilarityVector const& similarity_vector, MdLhs const& cur_node_lhs);
    void CollectRefinersForViolated(MdNode& cur_node, std::vector<MdRefiner>& found,
                                    MdLhs& cur_node_lhs, SimilarityVector const& similarity_vector,
                                    model::Index cur_node_index);

    bool IsUnsupported(MdLhs const& lhs) const;

    bool IsUnsupported(LhsSpecialization const& lhs_specialization) const;

    void UpdateMaxLevel(LhsSpecialization const& lhs_specialization);
    void AddNewMinimal(MdNode& cur_node, MdSpecialization const& md, model::Index cur_node_index);
    MdNode* TryGetNextNode(MdSpecialization const& md, GeneralizationHelper& helper,
                           model::Index cur_node_index, model::Index const next_node_index,
                           model::md::DecisionBoundary const next_lhs_bound);
    void AddIfMinimal(MdSpecialization const& md);

    static auto SetUnsupAction() noexcept {
        return [](SupportNode* node) { node->is_unsupported = true; };
    }

    // Generalization check, specialization (add if minimal)
    void MarkNewLhs(SupportNode& cur_node, MdLhs const& lhs, model::Index cur_node_index);
    void MarkUnsupported(MdLhs const& lhs);

    [[nodiscard]] std::optional<model::md::DecisionBoundary> SpecializeOneLhs(
            model::Index col_match_index, model::md::DecisionBoundary lhs_bound) const;
    void Specialize(MdLhs const& lhs, SimilarityVector const& specialize_past, Rhss const& rhss);
    void Specialize(MdLhs const& lhs, Rhss const& rhss);

    void GetAll(MdNode& cur_node, std::vector<MdLatticeNodeInfo>& collected, MdLhs& cur_node_lhs,
                model::Index this_node_index);

public:
    explicit MdLattice(std::size_t column_matches_size, SingleLevelFunc single_level_func,
                       std::vector<std::vector<model::md::DecisionBoundary>> const& lhs_bounds,
                       bool prune_nondisjoint);

    std::size_t GetColMatchNumber() const noexcept {
        return column_matches_size_;
    }

    [[nodiscard]] std::size_t GetMaxLevel() const noexcept {
        return max_level_;
    }

    std::vector<model::md::DecisionBoundary> GetRhsInterestingnessBounds(
            MdLhs const& lhs, std::vector<model::Index> const& indices) const;
    std::vector<MdVerificationMessenger> GetLevel(std::size_t level);
    std::vector<MdRefiner> CollectRefinersForViolated(SimilarityVector const& similarity_vector);
    std::vector<MdLatticeNodeInfo> GetAll();
};

}  // namespace algos::hymd::lattice

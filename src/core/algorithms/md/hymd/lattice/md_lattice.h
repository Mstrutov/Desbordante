#pragma once

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/column_match_info.h"
#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lattice/lattice_child_array.h"
#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/single_level_func.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/rhss.h"
#include "algorithms/md/hymd/similarity_vector.h"
#include "algorithms/md/hymd/utility/invalidated_rhss.h"
#include "model/index.h"

namespace algos::hymd::lattice {

class MdLattice {
private:
    struct MdNode {
        LatticeChildArray<MdNode> children;
        DecisionBoundaryVector rhs_bounds;

        MdNode(std::size_t attributes_num, std::size_t children_number)
            : children(children_number), rhs_bounds(attributes_num) {}

        explicit MdNode(DecisionBoundaryVector rhs)
            : children(rhs.size()), rhs_bounds(std::move(rhs)) {}
    };

    struct SupportNode {
        LatticeChildArray<SupportNode> children;
        bool is_unsupported = false;

        SupportNode(std::size_t children_number) : children(children_number) {}
    };

    class GeneralizationChecker;

    using MdBoundMap = BoundaryMap<MdNode>;
    using MdOptionalChild = std::optional<MdBoundMap>;
    using MdNodeChildren = LatticeChildArray<MdNode>;

    using SupportNodeChildren = LatticeChildArray<SupportNode>;

public:
    class MdRefiner {
        MdLattice* lattice_;
        SimilarityVector const* sim_;
        DecisionBoundaryVector lhs_;
        DecisionBoundaryVector* rhs_;
        utility::InvalidatedRhss invalidated_;

    public:
        MdRefiner(MdLattice* lattice, SimilarityVector const* sim, DecisionBoundaryVector lhs,
                  DecisionBoundaryVector* rhs, utility::InvalidatedRhss invalidated)
            : lattice_(lattice),
              sim_(sim),
              lhs_(std::move(lhs)),
              rhs_(rhs),
              invalidated_(std::move(invalidated)) {}

        void Refine();
    };

    class MdVerificationMessenger {
        MdLattice* lattice_;
        DecisionBoundaryVector lhs_;
        DecisionBoundaryVector* rhs_;

    public:
        MdVerificationMessenger(MdLattice* lattice, DecisionBoundaryVector lhs,
                                DecisionBoundaryVector* rhs)
            : lattice_(lattice), lhs_(std::move(lhs)), rhs_(rhs) {}

        DecisionBoundaryVector const& GetLhs() const {
            return lhs_;
        }

        DecisionBoundaryVector& GetRhs() {
            return *rhs_;
        }

        void MarkUnsupported();

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
    std::vector<ColumnMatchInfo> const* const column_matches_info_;
    bool const prune_nondisjoint_;

    bool HasLhsGeneralization(MdNode const& node, DecisionBoundaryVector const& lhs_bounds,
                              MdElement rhs, model::Index node_index,
                              model::Index start_index) const;

    void GetLevel(MdNode& cur_node, std::vector<MdVerificationMessenger>& collected,
                  DecisionBoundaryVector& cur_node_lhs_bounds, model::Index cur_node_index,
                  std::size_t level_left);

    [[nodiscard]] bool HasGeneralization(MdNode const& node,
                                         DecisionBoundaryVector const& lhs_bounds, MdElement rhs,
                                         model::Index cur_node_index) const;

    void RaiseInterestingnessBounds(
            MdNode const& cur_node, DecisionBoundaryVector const& lhs_bounds,
            std::vector<model::md::DecisionBoundary>& cur_interestingness_bounds,
            model::Index cur_node_index, std::vector<model::Index> const& indices) const;

    void TryAddRefiner(std::vector<MdRefiner>& found, DecisionBoundaryVector& rhs,
                       SimilarityVector const& similarity_vector,
                       DecisionBoundaryVector const& cur_node_lhs_bounds);
    void CollectRefinersForViolated(MdNode& cur_node, std::vector<MdRefiner>& found,
                                    DecisionBoundaryVector& cur_node_lhs_bounds,
                                    SimilarityVector const& similarity_vector,
                                    model::Index cur_node_index);

    void GetAll(MdNode& cur_node, std::vector<MdLatticeNodeInfo>& collected,
                DecisionBoundaryVector& cur_node_lhs_bounds, model::Index this_node_index);

    void AddNewMinimal(MdNode& cur_node, DecisionBoundaryVector const& lhs_bounds, MdElement rhs,
                       model::Index cur_node_index);

    void UpdateMaxLevel(DecisionBoundaryVector const& lhs_bounds);

    MdNode* ReturnNextNode(DecisionBoundaryVector const& lhs_bounds, GeneralizationChecker& checker,
                           model::Index cur_node_index, model::Index next_node_index);

    bool IsUnsupported(SupportNode const& cur_node, DecisionBoundaryVector const& lhs_bounds,
                       model::Index cur_node_index) const;
    void MarkNewLhs(SupportNode& cur_node, DecisionBoundaryVector const& lhs_bounds,
                    model::Index cur_node_index);

    [[nodiscard]] std::optional<model::md::DecisionBoundary> SpecializeOneLhs(
            model::Index col_match_index, model::md::DecisionBoundary lhs_bound) const;

    void AddIfMinimal(DecisionBoundaryVector const& lhs_bounds, MdElement rhs);
    bool IsUnsupported(DecisionBoundaryVector const& lhs_bounds) const;

    [[nodiscard]] bool HasGeneralization(DecisionBoundaryVector const& lhs_bounds,
                                         MdElement rhs) const;

    void MarkUnsupported(DecisionBoundaryVector const& lhs_bounds);

    void Specialize(DecisionBoundaryVector& lhs_bounds,
                    DecisionBoundaryVector const& specialize_past, Rhss const& rhss);

public:
    std::size_t GetColMatchNumber() const noexcept {
        return column_matches_size_;
    }

    [[nodiscard]] std::size_t GetMaxLevel() const noexcept {
        return max_level_;
    }

    std::vector<MdVerificationMessenger> GetLevel(std::size_t level);
    std::vector<model::md::DecisionBoundary> GetRhsInterestingnessBounds(
            DecisionBoundaryVector const& lhs_bounds,
            std::vector<model::Index> const& indices) const;
    std::vector<MdRefiner> CollectRefinersForViolated(SimilarityVector const& similarity_vector);
    std::vector<MdLatticeNodeInfo> GetAll();

    explicit MdLattice(std::size_t column_matches_size, SingleLevelFunc single_level_func,
                       std::vector<ColumnMatchInfo> const& column_matches_info,
                       bool prune_nondisjoint);
};

}  // namespace algos::hymd::lattice

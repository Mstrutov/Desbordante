#pragma once

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/node_base.h"
#include "algorithms/md/hymd/lattice/single_level_func.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/md_lhs.h"
#include "algorithms/md/hymd/rhss.h"
#include "algorithms/md/hymd/similarity_vector.h"
#include "algorithms/md/hymd/utility/invalidated_rhss.h"
#include "model/index.h"

namespace algos::hymd::lattice {

class MdLattice {
private:
    struct LhsSpecialization {
        MdLhs const& old_lhs;
        MdElement specialized;

        MdLhs const& ToUnspecialized() const {
            return old_lhs;
        }

        LhsSpecialization const& GetLhsSpecialization() const {
            return *this;
        }
    };

    struct Md {
        MdLhs const& lhs;
        MdElement rhs;
    };

    struct MdSpecialization {
        LhsSpecialization const& lhs_specialization;
        MdElement rhs;

        Md ToUnspecialized() const {
            return {lhs_specialization.ToUnspecialized(), rhs};
        }

        LhsSpecialization const& GetLhsSpecialization() const {
            return lhs_specialization;
        }
    };

    class MdNode : public NodeBase<MdNode> {
    public:
        using Specialization = MdSpecialization;
        using Unspecialized = Md;

        DecisionBoundaryVector rhs_bounds;

        static MdLhs const& GetLhs(Unspecialized const& md) {
            return md.lhs;
        }

        MdNode* AddOneUnchecked(model::Index child_array_index, model::md::DecisionBoundary bound) {
            return AddOneUncheckedBase(child_array_index, bound, rhs_bounds.size());
        }

        MdNode(std::size_t attributes_num, std::size_t children_number)
            : NodeBase<MdNode>(children_number), rhs_bounds(attributes_num) {}

        explicit MdNode(DecisionBoundaryVector rhs)
            : NodeBase<MdNode>(rhs.size()), rhs_bounds(std::move(rhs)) {}
    };

    class SupportNode : public NodeBase<SupportNode> {
    public:
        using Specialization = LhsSpecialization;
        using Unspecialized = MdLhs;

        bool is_unsupported = false;

        static MdLhs const& GetLhs(Unspecialized const& lhs) {
            return lhs;
        }

        SupportNode* AddOneUnchecked(model::Index child_array_index,
                                     model::md::DecisionBoundary bound) {
            return AddOneUncheckedBase(child_array_index, bound);
        }

        SupportNode(std::size_t children_number) : NodeBase<SupportNode>(children_number) {}
    };

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
        MdLhs lhs_;
        DecisionBoundaryVector* rhs_;
        utility::InvalidatedRhss invalidated_;

    public:
        MdRefiner(MdLattice* lattice, SimilarityVector const* sim, MdLhs lhs,
                  DecisionBoundaryVector* rhs, utility::InvalidatedRhss invalidated)
            : lattice_(lattice),
              sim_(sim),
              lhs_(std::move(lhs)),
              rhs_(rhs),
              invalidated_(std::move(invalidated)) {}

        MdLhs const& GetLhs() const {
            return lhs_;
        }

        void ZeroRhs() {
            rhs_->assign(lattice_->GetColMatchNumber(), 0.0);
        }

        void Refine();
    };

    class MdVerificationMessenger {
        MdLattice* lattice_;
        MdLhs lhs_;
        DecisionBoundaryVector* rhs_;

    public:
        MdVerificationMessenger(MdLattice* lattice, MdLhs lhs, DecisionBoundaryVector* rhs)
            : lattice_(lattice), lhs_(std::move(lhs)), rhs_(rhs) {}

        MdLhs const& GetLhs() const {
            return lhs_;
        }

        DecisionBoundaryVector& GetRhs() {
            return *rhs_;
        }

        void MarkUnsupported();

        void ZeroRhs() {
            rhs_->assign(lattice_->GetColMatchNumber(), 0.0);
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

    template <typename NodeType>
    bool HasChildGenSpec(NodeType const& node, model::Index node_index,
                         model::Index next_node_index, model::md::DecisionBoundary bound_limit,
                         auto const& md, auto gen_method, auto get_b_map_iter) const;
    template <typename NodeType>
    bool NodeHasLhsGeneralizationSpec(NodeType const& node,
                                      NodeType::Specialization const& specialization,
                                      model::Index node_index, model::Index start_index,
                                      auto spec_method, auto total_method) const;

    bool HasLhsGeneralizationSpec(MdNode const& node, MdSpecialization const& md,
                                  model::Index node_index, model::Index start_index) const;

    [[nodiscard]] bool HasGeneralizationSpec(MdNode const& node, MdSpecialization const& md,
                                             model::Index cur_node_index) const {
        return HasLhsGeneralizationSpec(node, md, cur_node_index, cur_node_index);
    }

    bool HasLhsGeneralizationTotal(MdNode const& node, Md const& md, model::Index node_index,
                                   model::Index start_index) const;

    [[nodiscard]] bool HasGeneralizationTotal(MdNode const& node, Md const& md,
                                              model::Index cur_node_index) const {
        if (node.rhs_bounds[md.rhs.index] >= md.rhs.decision_boundary) return true;
        return HasLhsGeneralizationTotal(node, md, cur_node_index, cur_node_index);
    }

    [[nodiscard]] bool HasGeneralization(Md const& md) const {
        return HasGeneralizationTotal(md_root_, md, 0);
    }

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

    bool IsUnsupportedTotal(SupportNode const& node, MdLhs const& lhs,
                            model::Index node_index) const;

    bool IsUnsupported(MdLhs const& lhs) const {
        return IsUnsupportedTotal(support_root_, lhs, 0);
    }

    bool IsUnsupportedSpec(SupportNode const& node, LhsSpecialization const& lhs_specialization,
                           model::Index node_index) const;

    bool IsUnsupported(LhsSpecialization const& lhs_specialization) const {
        return IsUnsupportedSpec(support_root_, lhs_specialization, 0);
    }

    void UpdateMaxLevel(LhsSpecialization const& lhs_specialization);
    void AddNewMinimal(MdNode& cur_node, MdSpecialization const& md, model::Index cur_node_index);
    MdNode* TryGetNextNode(MdSpecialization const& md, GeneralizationHelper& checker,
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

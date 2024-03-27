#pragma once

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lattice/lattice_child_array.h"
#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/single_level_func.h"
#include "algorithms/md/hymd/similarity_vector.h"
#include "model/index.h"

namespace algos::hymd::lattice {

class MdLattice {
private:
    struct Node {
        LatticeChildArray<Node> children;
        DecisionBoundaryVector rhs_bounds;

        Node(std::size_t attributes_num, std::size_t children_number)
            : children(children_number), rhs_bounds(attributes_num) {}

        explicit Node(DecisionBoundaryVector rhs)
            : children(rhs.size()), rhs_bounds(std::move(rhs)) {}
    };
    class GeneralizationChecker;

    using BoundMap = BoundaryMap<Node>;
    using OptionalChild = std::optional<BoundMap>;
    using NodeChildren = LatticeChildArray<Node>;

    std::size_t max_level_ = 0;
    Node root_;
    std::size_t const column_matches_size_;
    // Is there a way to define a level in such a way that one cannot use each decision boundary
    // independently to determine an MD's level but the lattice traversal algorithms still works?
    SingleLevelFunc const get_single_level_;

    bool HasLhsGeneralization(Node const& node, DecisionBoundaryVector const& lhs_bounds,
                              model::md::DecisionBoundary rhs_bound, model::Index rhs_index,
                              model::Index node_index, model::Index start_index) const;

    void GetLevel(Node& cur_node, std::vector<MdLatticeNodeInfo>& collected,
                  DecisionBoundaryVector& cur_node_lhs_bounds, model::Index cur_node_index,
                  std::size_t level_left);

    [[nodiscard]] bool HasGeneralization(Node const& node, DecisionBoundaryVector const& lhs_bounds,
                                         model::md::DecisionBoundary rhs_bound,
                                         model::Index rhs_index, model::Index cur_node_index) const;

    void RaiseInterestingnessBounds(
            Node const& cur_node, DecisionBoundaryVector const& lhs_bounds,
            std::vector<model::md::DecisionBoundary>& cur_interestingness_bounds,
            model::Index this_node_index, std::vector<model::Index> const& indices) const;

    void FindViolatedInternal(Node& cur_node, std::vector<MdLatticeNodeInfo>& found,
                              DecisionBoundaryVector& cur_node_lhs_bounds,
                              SimilarityVector const& similarity_vector,
                              model::Index this_node_index);

    void GetAll(Node& cur_node, std::vector<MdLatticeNodeInfo>& collected,
                DecisionBoundaryVector& cur_node_lhs_bounds, model::Index this_node_index);

    void AddNewMinimal(Node& cur_node, DecisionBoundaryVector const& lhs_bounds,
                       model::md::DecisionBoundary rhs_bound, model::Index rhs_index,
                       model::Index cur_node_index);

    void UpdateMaxLevel(DecisionBoundaryVector const& lhs_bounds);

    Node* ReturnNextNode(DecisionBoundaryVector const& lhs_bounds, GeneralizationChecker& checker,
                         model::Index cur_node_index, model::Index next_node_index);

public:
    std::size_t GetColMatchNumber() const noexcept {
        return column_matches_size_;
    }

    [[nodiscard]] bool HasGeneralization(DecisionBoundaryVector const& lhs_bounds,
                                         model::md::DecisionBoundary rhs_bound,
                                         model::Index rhs_index) const;

    [[nodiscard]] std::size_t GetMaxLevel() const noexcept {
        return max_level_;
    }

    std::vector<MdLatticeNodeInfo> GetLevel(std::size_t level);
    std::vector<model::md::DecisionBoundary> GetRhsInterestingnessBounds(
            DecisionBoundaryVector const& lhs_bounds,
            std::vector<model::Index> const& indices) const;
    void AddIfMinimal(DecisionBoundaryVector const& lhs_bounds,
                      model::md::DecisionBoundary rhs_bound, model::Index rhs_index);
    std::vector<MdLatticeNodeInfo> FindViolated(SimilarityVector const& similarity_vector);
    std::vector<MdLatticeNodeInfo> GetAll();

    explicit MdLattice(std::size_t column_matches_size, SingleLevelFunc single_level_func);
};

}  // namespace algos::hymd::lattice

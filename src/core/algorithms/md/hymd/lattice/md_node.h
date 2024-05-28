#pragma once

#include "algorithms/md/hymd/lattice/md.h"
#include "algorithms/md/hymd/lattice/md_specialization.h"
#include "algorithms/md/hymd/lattice/node_base.h"
#include "algorithms/md/hymd/lattice/rhs.h"

namespace algos::hymd::lattice {
class MdNode : public NodeBase<MdNode> {
public:
    using Specialization = MdSpecialization;

    Rhs rhs;

    static MdLhs const& GetLhs(Specialization::Unspecialized const& md) noexcept {
        return md.lhs;
    }

    bool ContainsGeneralizationOf(Specialization::Unspecialized const& md) const noexcept {
        return rhs[md.rhs.index] >= md.rhs.decision_boundary;
    }

    MdNode* AddOneUnchecked(model::Index child_array_index, model::md::DecisionBoundary bound) {
        return AddOneUncheckedBase(child_array_index, bound, rhs.size());
    }

    MdNode(std::size_t attributes_num, std::size_t children_number)
        : NodeBase<MdNode>(children_number), rhs(attributes_num) {}

    explicit MdNode(Rhs rhs) : NodeBase<MdNode>(rhs.size()), rhs(std::move(rhs)) {}
};
}  // namespace algos::hymd::lattice

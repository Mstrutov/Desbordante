//
// Created by alexandrsmirn
//

#pragma once

#include "Vertical.h"
#include "DFD/LatticeObservations/LatticeObservations.h"

class NonDependenciesMap : public std::unordered_map<Vertical, std::unordered_set<Vertical>> {
    using vertical_set = std::unordered_set<Vertical>;
public:

    explicit NonDependenciesMap(RelationalSchema const* schema);
    NonDependenciesMap() = default;

    std::unordered_set<Vertical> getPrunedSupersets(std::unordered_set<Vertical> const& supersets) const;
    void addNewNonDependency(Vertical const& node);
    bool canBePruned(Vertical const& node) const;
    void rebalance();
};


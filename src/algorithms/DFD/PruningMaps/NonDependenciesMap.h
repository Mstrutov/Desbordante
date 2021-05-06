//
// Created by alexandrsmirn
//

#pragma once

#include "Vertical.h"
#include "DFD/LatticeObservations/LatticeObservations.h"
#include "CustomComparator.h"

class NonDependenciesMap : public std::unordered_map<Vertical, std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>> {
    using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;
public:

    explicit NonDependenciesMap(shared_ptr<RelationalSchema> schema);
    NonDependenciesMap() = default;

    vertical_set getPrunedSupersets(vertical_set const& supersets) const;
    void addNewNonDependency(shared_ptr<Vertical> const& node);
    bool canBePruned(Vertical const& node) const;
    void rebalance();
};


//
// Created by alexandrsmirn
//

#pragma once

#include <unordered_map>
#include <unordered_set>

#include "Vertical.h"
#include "DFD/LatticeObservations/LatticeObservations.h"

class DependenciesMap : public std::unordered_map<Vertical, std::unordered_set<Vertical>> { //TODO указателями или просто объектами?

public:
    using vertical_set = std::unordered_set<Vertical>;
    explicit DependenciesMap(RelationalSchema const* schema);
    DependenciesMap() = default;

    void addNewDependency(Vertical const& nodeToAdd);
    //vector<shared_ptr<Vertical>> getUncheckedSubsets(shared_ptr<Vertical> node, LatticeObservations & observations);
    std::unordered_set<Vertical> getPrunedSubsets(std::unordered_set<Vertical> const& subsets) const;
    bool canBePruned(Vertical const& node) const;
    void rebalance();
};

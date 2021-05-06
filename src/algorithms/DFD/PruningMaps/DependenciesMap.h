//
// Created by alexandrsmirn
//

#pragma once

#include <unordered_map>
#include <unordered_set>

#include "Vertical.h"
#include "DFD/LatticeObservations/LatticeObservations.h"

class DependenciesMap : public std::unordered_map<Vertical, std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>> { //TODO указателями или просто объектами?

public:
    using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;
    explicit DependenciesMap(shared_ptr<RelationalSchema> schema);
    DependenciesMap() = default;

    void addNewDependency(shared_ptr<Vertical> nodeToAdd);
    //vector<shared_ptr<Vertical>> getUncheckedSubsets(shared_ptr<Vertical> node, LatticeObservations & observations);
    vertical_set getPrunedSubsets(vertical_set const& subsets) const;
    bool canBePruned(Vertical const& node) const;
    void rebalance();
};

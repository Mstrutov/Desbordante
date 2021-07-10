//
// Created by alexandrsmirn
//

#pragma once

#include "Vertical.h"
#include "PruningMap.h"

class DependenciesMap : public PruningMap { //TODO указателями или просто объектами?
public:
    explicit DependenciesMap(RelationalSchema const* schema);
    DependenciesMap() = default;

    void addNewDependency(Vertical const& nodeToAdd);
    std::unordered_set<Vertical> getPrunedSubsets(std::unordered_set<Vertical> const& subsets) const;
    bool canBePruned(Vertical const& node) const;
};

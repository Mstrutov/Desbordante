#pragma once

#include "Vertical.h"
#include "PruningMap.h"

class DependenciesMap : public PruningMap {
public:
    explicit DependenciesMap(RelationalSchema const* schema);
    DependenciesMap() = default;

    std::unordered_set<Vertical> getPrunedSubsets(std::unordered_set<Vertical> const& subsets) const;
    void addNewDependency(Vertical const& nodeToAdd);
    bool canBePruned(Vertical const& node) const;
};

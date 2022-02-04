#pragma once

#include "Vertical.h"
#include "PruningMap.h"

class DependenciesMap : public PruningMap {
public:
    explicit DependenciesMap(RelationalSchema const* schema);
    DependenciesMap() = default;

    std::unordered_set<Vertical> GetPrunedSubsets(std::unordered_set<Vertical> const& subsets) const;
    void AddNewDependency(Vertical const& node_to_add);
    bool CanBePruned(Vertical const& node) const;
};

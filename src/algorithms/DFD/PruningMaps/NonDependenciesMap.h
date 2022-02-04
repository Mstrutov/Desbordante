#pragma once

#include "Vertical.h"
#include "PruningMap.h"

class NonDependenciesMap : public PruningMap {
public:
    explicit NonDependenciesMap(RelationalSchema const* schema);
    NonDependenciesMap() = default;

    std::unordered_set<Vertical> GetPrunedSupersets(std::unordered_set<Vertical> const& supersets) const;
    void AddNewNonDependency(Vertical const& node_to_add);
    bool CanBePruned(Vertical const& node) const;
};

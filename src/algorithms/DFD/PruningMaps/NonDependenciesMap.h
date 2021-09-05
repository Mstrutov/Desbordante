#pragma once

#include "Vertical.h"
#include "PruningMap.h"

class NonDependenciesMap : public PruningMap {
public:
    explicit NonDependenciesMap(RelationalSchema const* schema);
    NonDependenciesMap() = default;

    std::unordered_set<Vertical> getPrunedSupersets(std::unordered_set<Vertical> const& supersets) const;
    void addNewNonDependency(Vertical const& node);
    bool canBePruned(Vertical const& node) const;
};

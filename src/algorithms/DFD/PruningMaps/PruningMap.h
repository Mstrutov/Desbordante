#pragma once

#include <unordered_set>
#include <unordered_map>

#include "Vertical.h"
#include "DFD/LatticeObservations/LatticeObservations.h"

class PruningMap : public std::unordered_map<Vertical, std::unordered_set<Vertical>> {
public:
    PruningMap(RelationalSchema const* schema);
    PruningMap() = default;

    void rebalance();
    void rebalanceGroup(Vertical const& key);
};

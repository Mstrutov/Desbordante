//
// Created by alexandrsmirn
//

#pragma once

#include <unordered_map>
#include <unordered_set>

#include "Vertical.h"
#include "DFD/LatticeObservations/LatticeObservations.h"

/*namespace std {
    template<>
    struct hash<Vertical> {
        size_t operator()(Vertical const &k) const {
            return k.getColumnIndices().to_ulong();
        }
    };
}*/

/*namespace std {
    template<>
    struct hash<shared_ptr<Vertical>> {
        size_t operator()(shared_ptr<Vertical> const &k) const {
            return k->getColumnIndices().to_ulong();
        }
    };
}*/

class DependenciesMap : public std::unordered_map<Vertical, std::unordered_set<shared_ptr<Vertical>>> { //TODO указателями или просто объектами?

public:
    explicit DependenciesMap(shared_ptr<RelationalSchema> schema);
    DependenciesMap() = default;

    void addNewDependency(shared_ptr<Vertical> node);
    vector<shared_ptr<Vertical>> getUncheckedSubsets(shared_ptr<Vertical> node, LatticeObservations const& observations);
    bool canBePruned(Vertical const& node);
};

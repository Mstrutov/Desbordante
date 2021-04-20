//
// Created by alexandrsmirn
//

#pragma once

#include "Vertical.h"
#include "DFD/LatticeObservations/LatticeObservations.h"

namespace std {
    template<>
    struct hash<shared_ptr<Vertical>> {
        size_t operator()(shared_ptr<Vertical> const &k) const {
            return k->getColumnIndices().to_ulong();
        }
    };
}

class NonDependenciesMap : public std::unordered_map<Vertical, std::unordered_set<shared_ptr<Vertical>>>{
public:
    explicit NonDependenciesMap(shared_ptr<RelationalSchema> schema);
    NonDependenciesMap() = default;

    void addNewNonDependency(shared_ptr<Vertical> node);
    vector<shared_ptr<Vertical>> getUncheckedSupersets(shared_ptr<Vertical> node, LatticeObservations const& observations);
    bool canBePruned(Vertical const& node);
};


#include "NonDependenciesMap.h"

NonDependenciesMap::NonDependenciesMap(RelationalSchema const* schema)
    : PruningMap(schema) {}

std::unordered_set<Vertical> NonDependenciesMap::GetPrunedSupersets(std::unordered_set<Vertical> const& supersets) const {
    std::unordered_set<Vertical> pruned_supersets;
    for (auto const& node : supersets) {
        if (CanBePruned(node)) {
            pruned_supersets.insert(node);
        }
    }
    return pruned_supersets;
}

bool NonDependenciesMap::CanBePruned(const Vertical& node) const {
    for (auto const& map_row : *this) {
        Vertical const& key = map_row.first;
        if (node.Contains(key)) {
            for (Vertical const& non_dependency : map_row.second) {
                if (non_dependency.Contains(node)) {
                    return true;
                }
            }
        }
    }
    return false;
}

void NonDependenciesMap::AddNewNonDependency(Vertical const& node_to_add) {
    for (auto& map_row : *this) {
        Vertical const& key = map_row.first;

        if (node_to_add.Contains(key)) {
            auto& non_deps_for_key = map_row.second;
            bool has_superset_entry = false;

            for (auto iter = non_deps_for_key.begin(); iter != non_deps_for_key.end();) {
                //if verticals are the same, then contains == true
                Vertical const& non_dep = *iter;
                if (non_dep.Contains(node_to_add)) {
                    has_superset_entry = true;
                    break;
                } else if (node_to_add.Contains(non_dep)) {
                    iter = non_deps_for_key.erase(iter);
                } else {
                    iter++;
                }
            }

            if (!has_superset_entry) {
                non_deps_for_key.insert(node_to_add);
            }
        }
    }
    Rebalance();
}

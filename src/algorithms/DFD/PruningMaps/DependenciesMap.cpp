#include "DependenciesMap.h"

DependenciesMap::DependenciesMap(RelationalSchema const* schema)
    : PruningMap(schema) {}

std::unordered_set<Vertical> DependenciesMap::GetPrunedSubsets(std::unordered_set<Vertical> const& subsets) const {
    std::unordered_set<Vertical> pruned_subsets;
    for (auto const& node : subsets) {
        if (CanBePruned(node)) {
            pruned_subsets.insert(node);
        }
    }
    return pruned_subsets;
}

void DependenciesMap::AddNewDependency(Vertical const& node_to_add) {
    for (auto& map_row : *this) {
        Vertical const& key = map_row.first;

        if (node_to_add.Contains(key)) {
            auto& deps_for_key = map_row.second;
            bool has_subset_entry = false;

            for (auto iter = deps_for_key.begin(); iter != deps_for_key.end();) {
                //if verticals are the same, then contains == true
                Vertical const& dep = *iter;
                if (node_to_add.Contains(dep)) {
                    has_subset_entry = true;
                    break;
                } else if (dep.Contains(node_to_add)) {
                    iter = deps_for_key.erase(iter);
                } else {
                    iter++;
                }
            }

            if (!has_subset_entry) {
                deps_for_key.insert(node_to_add);
            }
        }
    }
    Rebalance();
}

bool DependenciesMap::CanBePruned(Vertical const& node) const {
    for (auto const& map_row : *this) {
        Vertical const& key = map_row.first;
        if (node.Contains(key)) {
            for (Vertical const& dependency : map_row.second) {
                if (node.Contains(dependency)) {
                    return true;
                }
            }
        }
    }
    return false;
}

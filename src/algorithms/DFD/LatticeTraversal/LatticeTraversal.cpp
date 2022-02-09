#include "LatticeTraversal.h"

#include <random>

#include "PositionListIndex.h"

LatticeTraversal::LatticeTraversal(const Column* const rhs,
                                   const ColumnLayoutRelationData* const relation,
                                   const std::vector<Vertical>& unique_verticals,
                                   PartitionStorage* const partition_storage)
    : rhs_(rhs),
      dependencies_map_(relation->GetSchema()),
      non_dependencies_map_(relation->GetSchema()),
      column_order_(relation),
      unique_columns_(unique_verticals),
      relation_(relation),
      partition_storage_(partition_storage),
      gen_(rd_()) {}

std::unordered_set<Vertical> LatticeTraversal::FindLHSs() {
    RelationalSchema const* const schema = relation_->GetSchema();

    //processing of found unique columns
    for (auto const& lhs : unique_columns_) {
        if (!lhs.Contains(*rhs_)) {
            observations_[lhs] = NodeCategory::kMinimalDependency;
            dependencies_map_.AddNewDependency(lhs);
            minimal_deps_.insert(lhs);
        }
    }

    std::stack<Vertical> seeds;

    /* Temporary fix. I think `GetOrderHighDistinctCount` should return vector of
     * unsigned integers since `order` sould be something non-negative.
     */
    for (unsigned partition_index :
         column_order_.GetOrderHighDistinctCount(Vertical(*rhs_).Invert())) {
        if (partition_index != rhs_->GetIndex()) {
            seeds.push(Vertical(*schema->GetColumn(partition_index)));
        }
    }

    do {
        while (!seeds.empty()) {
            Vertical node;
            if (!seeds.empty()) {
                node = std::move(seeds.top());
                seeds.pop();
            } else {
                node = *schema->empty_vertical_;
            }

            do {
                auto const node_observation_iter = observations_.find(node);

                if (node_observation_iter != observations_.end()) {
                    NodeCategory& node_category = node_observation_iter->second;

                    if (node_category == NodeCategory::kCandidateMinimalDependency) {
                        node_category = observations_.UpdateDependencyCategory(node);
                        if (node_category == NodeCategory::kMinimalDependency) {
                            minimal_deps_.insert(node);
                        }
                    } else if (node_category == NodeCategory::kCandidateMaximalNonDependency) {
                        node_category =
                            observations_.UpdateNonDependencyCategory(node, rhs_->GetIndex());
                        if (node_category == NodeCategory::kMaximalNonDependency) {
                            maximal_non_deps_.insert(node);
                        }
                    }
                } else if (!InferCategory(node, rhs_->GetIndex())) {
                    //if we were not able to infer category, we calculate the partitions
                    auto node_pli = partition_storage_->GetOrCreateFor(node);
                    auto node_pli_pointer =
                        std::holds_alternative<util::PositionListIndex*>(node_pli)
                            ? std::get<util::PositionListIndex*>(node_pli)
                            : std::get<std::unique_ptr<util::PositionListIndex>>(node_pli).get();
                    auto intersected_pli = partition_storage_->GetOrCreateFor(node.Union(*rhs_));
                    auto intersected_pli_pointer =
                        std::holds_alternative<util::PositionListIndex*>(intersected_pli)
                            ? std::get<util::PositionListIndex*>(intersected_pli)
                            : std::get<std::unique_ptr<util::PositionListIndex>>(intersected_pli)
                                  .get();

                    if (node_pli_pointer->GetNepAsLong() ==
                        intersected_pli_pointer->GetNepAsLong()) {
                        observations_.UpdateDependencyCategory(node);
                        if (observations_[node] == NodeCategory::kMinimalDependency) {
                            minimal_deps_.insert(node);
                        }
                        dependencies_map_.AddNewDependency(node);
                    } else {
                        observations_.UpdateNonDependencyCategory(node, rhs_->GetIndex());
                        if (observations_[node] == NodeCategory::kMaximalNonDependency) {
                            maximal_non_deps_.insert(node);
                        }
                        non_dependencies_map_.AddNewNonDependency(node);
                    }
                }

                node = PickNextNode(node, rhs_->GetIndex());
            } while (node != *node.GetSchema()->empty_vertical_);
        }
        seeds = GenerateNextSeeds(rhs_);
    } while (!seeds.empty());

    return minimal_deps_;
}

bool LatticeTraversal::InferCategory(Vertical const& node, unsigned int rhs_index) {
    if (non_dependencies_map_.CanBePruned(node)) {
        observations_.UpdateNonDependencyCategory(node, rhs_index);
        non_dependencies_map_.AddNewNonDependency(node);
        if (observations_[node] == NodeCategory::kMinimalDependency) {
            minimal_deps_.insert(node);
        }
        return true;
    } else if (dependencies_map_.CanBePruned(node)) {
        observations_.UpdateDependencyCategory(node);
        dependencies_map_.AddNewDependency(node);
        if (observations_[node] == NodeCategory::kMaximalNonDependency) {
            maximal_non_deps_.insert(node);
        }
        return true;
    }

    return false;
}

Vertical const& LatticeTraversal::TakeRandom(std::unordered_set<Vertical>& node_set) {
    std::uniform_int_distribution<> dis(0, std::distance(node_set.begin(), node_set.end()) - 1);
    auto iterator = node_set.begin();
    std::advance(iterator, dis(this->gen_));
    Vertical const& node = *iterator;
    return node;
}

Vertical LatticeTraversal::PickNextNode(Vertical const& node, unsigned int rhs_index) {
    auto node_iter = observations_.find(node);

    if (node_iter != observations_.end()) {
        if (node_iter->second == NodeCategory::kCandidateMinimalDependency) {
            auto unchecked_subsets = observations_.GetUncheckedSubsets(node, column_order_);
            auto pruned_non_dep_subsets =
                non_dependencies_map_.GetPrunedSupersets(unchecked_subsets);
            for (auto const& pruned_subset : pruned_non_dep_subsets) {
                observations_[pruned_subset] = NodeCategory::kNonDependency;
            }
            SubstractSets(unchecked_subsets, pruned_non_dep_subsets);

            if (unchecked_subsets.empty() && pruned_non_dep_subsets.empty()) {
                minimal_deps_.insert(node);
                observations_[node] = NodeCategory::kMinimalDependency;
            } else if (!unchecked_subsets.empty()) {
                auto const& next_node = TakeRandom(unchecked_subsets);
                //auto const& next_node = *unchecked_subsets.begin();
                trace_.push(node);
                return next_node;
            }
        } else if (node_iter->second == NodeCategory::kCandidateMaximalNonDependency) {
            auto unchecked_supersets =
                observations_.GetUncheckedSupersets(node, rhs_index, column_order_);
            auto pruned_non_dep_supersets =
                non_dependencies_map_.GetPrunedSupersets(unchecked_supersets);
            auto pruned_dep_supersets = dependencies_map_.GetPrunedSubsets(unchecked_supersets);

            for (auto const& pruned_superset : pruned_non_dep_supersets) {
                observations_[pruned_superset] = NodeCategory::kNonDependency;
            }
            for (auto const& pruned_superset : pruned_dep_supersets) {
                observations_[pruned_superset] = NodeCategory::kDependency;
            }

            SubstractSets(unchecked_supersets, pruned_dep_supersets);
            SubstractSets(unchecked_supersets, pruned_non_dep_supersets);

            if (unchecked_supersets.empty() && pruned_non_dep_supersets.empty()) {
                maximal_non_deps_.insert(node);
                observations_[node] = NodeCategory::kMaximalNonDependency;
            } else if (!unchecked_supersets.empty()) {
                auto const& next_node = TakeRandom(unchecked_supersets);
                //auto const& next_node = *unchecked_supersets.begin();
                trace_.push(node);
                return next_node;
            }
        }
    }

    Vertical next_node = *(node.GetSchema()->empty_vertical_);
    if (!trace_.empty()) {
        next_node = trace_.top();
        trace_.pop();
    }
    return next_node;
}

std::stack<Vertical> LatticeTraversal::GenerateNextSeeds(Column const* const current_rhs) {
    std::unordered_set<Vertical> seeds;
    std::unordered_set<Vertical> new_seeds;

    for (auto const& non_dep : maximal_non_deps_) {
        auto complement_indices = non_dep.GetColumnIndicesRef();
        complement_indices[current_rhs->GetIndex()] = true;
        complement_indices.flip();

        if (seeds.empty()) {
            boost::dynamic_bitset<> single_column_bitset(relation_->GetNumColumns(), 0);
            single_column_bitset.reset();

            for (size_t column_index = complement_indices.find_first();
                 column_index < complement_indices.size();
                 column_index = complement_indices.find_next(column_index)) {
                single_column_bitset[column_index] = true;
                seeds.emplace(relation_->GetSchema(), single_column_bitset);
                single_column_bitset[column_index] = false;
            }
        } else {
            for (auto const& dependency : seeds) {
                auto new_combination = dependency.GetColumnIndicesRef();

                for (size_t column_index = complement_indices.find_first();
                     column_index < complement_indices.size();
                     column_index = complement_indices.find_next(column_index)) {
                    new_combination[column_index] = true;
                    new_seeds.emplace(relation_->GetSchema(), new_combination);
                    new_combination[column_index] = dependency.GetColumnIndicesRef()[column_index];
                }
            }

            std::list<Vertical> minimized_new_seeds = Minimize(new_seeds);
            seeds.clear();
            for (auto& new_seed : minimized_new_seeds) {
                seeds.insert(std::move(new_seed));
            }
            new_seeds.clear();
        }
    }

    for (auto seed_iter = seeds.begin(); seed_iter != seeds.end();) {
        if (minimal_deps_.find(*seed_iter) != minimal_deps_.end()) {
            seed_iter = seeds.erase(seed_iter);
        } else {
            seed_iter++;
        }
    }

    std::stack<Vertical> remaining_seeds;

    for (auto const& new_seed : seeds) {
        remaining_seeds.push(new_seed);
    }

    return remaining_seeds;
}

std::list<Vertical> LatticeTraversal::Minimize(
    std::unordered_set<Vertical> const& node_list) const {
    unsigned int max_cardinality = 0;
    std::unordered_map<unsigned int, std::list<Vertical const*>> seeds_by_size(
        node_list.size() / relation_->GetNumColumns());

    for (auto const& seed : node_list) {
        unsigned int const cardinality_of_seed = seed.GetArity();
        max_cardinality = std::max(max_cardinality, cardinality_of_seed);
        if (seeds_by_size.find(cardinality_of_seed) == seeds_by_size.end()) {
            seeds_by_size[cardinality_of_seed] = std::list<Vertical const*>();
        }
        seeds_by_size[cardinality_of_seed].push_back(&seed);
    }

    for (unsigned int lower_bound = 1; lower_bound < max_cardinality; ++lower_bound) {
        if (seeds_by_size.find(lower_bound) != seeds_by_size.end()) {
            auto const& lower_bound_seeds = seeds_by_size.find(lower_bound)->second;
            for (unsigned int upper_bound = max_cardinality; upper_bound > lower_bound;
                 --upper_bound) {
                if (seeds_by_size.find(upper_bound) != seeds_by_size.end()) {
                    auto& upper_bound_seeds = seeds_by_size.find(upper_bound)->second;
                    for (auto const& lower_seed : lower_bound_seeds) {
                        for (auto upper_it = upper_bound_seeds.begin();
                             upper_it != upper_bound_seeds.end();) {
                            if ((*upper_it)->Contains(*lower_seed)) {
                                upper_it = upper_bound_seeds.erase(upper_it);
                            } else {
                                ++upper_it;
                            }
                        }
                    }
                }
            }
        }
    }

    std::list<Vertical> new_seeds;
    for (auto& seed_list : seeds_by_size) {
        for (auto& seed : seed_list.second) {
            new_seeds.push_back(*seed);
        }
    }
    return new_seeds;
}

void LatticeTraversal::SubstractSets(std::unordered_set<Vertical>& set,
                                     std::unordered_set<Vertical> const& set_to_substract) {
    for (const auto& node_to_delete : set_to_substract) {
        auto found_element_iter = set.find(node_to_delete);
        if (found_element_iter != set.end()) {
            set.erase(found_element_iter);
        }
    }
}

#include "order.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <utility>

#include <easylogging++.h>

#include "config/names_and_descriptions.h"
#include "config/tabular_data/input_table/option.h"
#include "dependency_checker.h"
#include "list_lattice.h"
#include "model/table/tuple_index.h"
#include "model/types/types.h"
#include "order_utility.h"

namespace algos::order {

Order::Order() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName()});
}

void Order::RegisterOptions() {
    using namespace config::names;
    using namespace config::descriptions;
    using config::Option;

    RegisterOption(config::kTableOpt(&input_table_));
}

void Order::LoadDataInternal() {
    typed_relation_ = model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, false);
}

void Order::ResetState() {}

void Order::PruneSingleEqClassPartitions() {
    for (auto& [attr, partition] : sorted_partitions_) {
        if (partition.Size() == 1) {
            for (AttributeList other_attr : single_attributes_) {
                if (other_attr != attr) {
                    valid_[attr].insert(other_attr);
                }
            }
            single_attributes_.erase(
                    std::find(single_attributes_.begin(), single_attributes_.end(), attr));
        }
    }
}

void Order::CreateSingleColumnSortedPartitions() {
    std::vector<model::TypedColumnData> const& data = typed_relation_->GetColumnData();
    std::unordered_set<model::TupleIndex> null_rows = GetNullIndices(data);
    for (unsigned int i = 0; i < data.size(); ++i) {
        if (!model::Type::IsOrdered(data[i].GetTypeId())) {
            continue;
        }
        single_attributes_.push_back({i});
        std::vector<IndexedByteData> indexed_byte_data = GetIndexedByteData(data[i], null_rows);
        std::unique_ptr<model::Type> type = model::CreateType(data[i].GetTypeId(), true);
        std::unique_ptr<model::MixedType> mixed_type =
                model::CreateSpecificType<model::MixedType>(model::TypeId::kMixed, true);
        auto less = [&type, &mixed_type](IndexedByteData const& l, IndexedByteData const& r) {
            if (type->GetTypeId() == +(model::TypeId::kMixed)) {
                return mixed_type->CompareAsStrings(l.data, r.data) == model::CompareResult::kLess;
            }
            return type->Compare(l.data, r.data) == model::CompareResult::kLess;
        };
        auto equal = [&type, &mixed_type](IndexedByteData const& l, IndexedByteData const& r) {
            if (type->GetTypeId() == +(model::TypeId::kMixed)) {
                return mixed_type->CompareAsStrings(l.data, r.data) == model::CompareResult::kEqual;
            }
            return type->Compare(l.data, r.data) == model::CompareResult::kEqual;
        };
        std::sort(indexed_byte_data.begin(), indexed_byte_data.end(), less);
        SortedPartition::EquivalenceClasses equivalence_classes;
        equivalence_classes.reserve(typed_relation_->GetNumRows());
        equivalence_classes.push_back({indexed_byte_data.front().index});
        for (size_t k = 1; k < indexed_byte_data.size(); ++k) {
            if (equal(indexed_byte_data[k - 1], indexed_byte_data[k])) {
                equivalence_classes.back().insert(indexed_byte_data[k].index);
            } else {
                equivalence_classes.push_back({indexed_byte_data[k].index});
            }
        }
        equivalence_classes.shrink_to_fit();
        sorted_partitions_.emplace(
                AttributeList{i},
                SortedPartition(std::move(equivalence_classes), typed_relation_->GetNumRows()));
    }
    PruneSingleEqClassPartitions();
}

void Order::CreateSortedPartitionsFromSingletons(AttributeList const& attr_list) {
    if (sorted_partitions_.find(attr_list) != sorted_partitions_.end()) {
        return;
    }
    SortedPartition res = sorted_partitions_.at({attr_list[0]});
    for (size_t i = 1; i < attr_list.size(); ++i) {
        res.Intersect(sorted_partitions_.at({attr_list[i]}));
    }
    sorted_partitions_.emplace(attr_list, res);
}

bool Order::HasValidPrefix(AttributeList const& lhs, AttributeList const& rhs) const {
    bool prefix_valid = false;
    for (AttributeList const& rhs_prefix : GetPrefixes(rhs)) {
        if (InUnorderedMap(valid_, lhs, rhs_prefix)) {
            prefix_valid = true;
            break;
        }
    }
    return prefix_valid;
}

ValidityType Order::CheckCandidateValidity(AttributeList const& lhs, AttributeList const& rhs) {
    bool is_merge_immediately = false;
    for (AttributeList const& lhs_prefix : GetPrefixes(lhs)) {
        if (InUnorderedMap(merge_invalidated_, lhs_prefix, rhs)) {
            is_merge_immediately = true;
            break;
        }
    }
    ValidityType candidate_validity = +ValidityType::merge;
    if (!is_merge_immediately) {
        CreateSortedPartitionsFromSingletons(lhs);
        if (sorted_partitions_[lhs].Size() == 1) {
            candidate_validity = +ValidityType::valid;
            candidate_sets_[lhs].erase(rhs);
        } else {
            CreateSortedPartitionsFromSingletons(rhs);
            candidate_validity = CheckForSwap(sorted_partitions_[lhs], sorted_partitions_[rhs]);
        }
    }
    return candidate_validity;
}

void Order::ComputeDependencies(ListLattice::LatticeLevel const& lattice_level) {
    if (lattice_->GetLevelNumber() < 2) {
        return;
    }
    UpdateCandidateSets();
    for (Node const& node : lattice_level) {
        CandidatePairs candidate_pairs = lattice_->ObtainCandidates(node);
        for (auto const& [lhs, rhs] : candidate_pairs) {
            if (!InUnorderedMap(candidate_sets_, lhs, rhs)) {
                continue;
            }
            if (HasValidPrefix(lhs, rhs)) {
                continue;
            }
            ValidityType candidate_validity = CheckCandidateValidity(lhs, rhs);
            if (candidate_validity == +ValidityType::valid) {
                bool non_minimal_by_merge = false;
                for (AttributeList const& merge_lhs : GetPrefixes(lhs)) {
                    if (InUnorderedMap(merge_invalidated_, merge_lhs, rhs)) {
                        non_minimal_by_merge = true;
                        break;
                    }
                }
                if (non_minimal_by_merge) {
                    continue;
                }
                if (valid_.find(lhs) == valid_.end()) {
                    valid_[lhs] = {};
                }
                valid_[lhs].insert(rhs);
                bool lhs_unique = typed_relation_->GetNumRows() == sorted_partitions_[lhs].Size();
                if (lhs_unique) {
                    candidate_sets_[lhs].erase(rhs);
                }
            } else if (candidate_validity == +ValidityType::swap) {
                candidate_sets_[lhs].erase(rhs);
            } else if (candidate_validity == +ValidityType::merge) {
                if (merge_invalidated_.find(lhs) == merge_invalidated_.end()) {
                    merge_invalidated_[lhs] = {};
                }
                merge_invalidated_[lhs].insert(rhs);
            }
        }
    }
    MergePrune();
}

std::vector<AttributeList> Order::Extend(AttributeList const& lhs, AttributeList const& rhs) const {
    std::vector<AttributeList> extended_rhss;
    for (auto const& single_attribute : single_attributes_) {
        if (AreDisjoint(single_attribute, lhs) && AreDisjoint(single_attribute, rhs)) {
            AttributeList extended(rhs);
            extended.push_back(single_attribute[0]);
            extended_rhss.push_back(extended);
        }
    }
    return extended_rhss;
}

bool Order::IsMinimal(AttributeList const& a) const {
    for (auto const& [lhs, rhs_list] : valid_) {
        for (AttributeList const& rhs : rhs_list) {
            auto it_rhs = std::search(a.begin(), a.end(), rhs.begin(), rhs.end());
            if (it_rhs == a.end()) {
                continue;
            }
            if (std::search(it_rhs + rhs.size(), a.end(), lhs.begin(), lhs.end()) != a.end()) {
                return false;
            }
            auto it_lhs = std::search(a.begin(), it_rhs, lhs.begin(), lhs.end());
            if (it_lhs + lhs.size() == it_rhs) {
                return false;
            }
        }
    }
    return true;
}

void Order::UpdateCandidateSets() {
    unsigned int level_num = lattice_->GetLevelNumber();
    if (level_num < 3) {
        return;
    }
    CandidateSets next_candidates;
    for (auto const& [lhs, rhs_list] : candidate_sets_) {
        next_candidates[lhs] = {};
        if (lhs.size() != level_num - 1) {
            for (AttributeList const& rhs : rhs_list) {
                if (InUnorderedMap(valid_, lhs, rhs)) {
                    continue;
                }
                std::vector<AttributeList> extended_rhss = Extend(lhs, rhs);
                for (AttributeList const& extended : extended_rhss) {
                    if (lhs.size() > 1) {
                        AttributeList lhs_max_prefix = MaxPrefix(lhs);
                        std::vector<AttributeList> extended_prefixes = GetPrefixes(extended);
                        auto prefix_is_valid =
                                [this, &lhs_max_prefix](AttributeList const& extended_prefix) {
                                    return InUnorderedMap(valid_, lhs_max_prefix, extended_prefix);
                                };
                        if ((std::find_if(extended_prefixes.begin(), extended_prefixes.end(),
                                          prefix_is_valid) == extended_prefixes.end()) &&
                            !InUnorderedMap(candidate_sets_, lhs_max_prefix, extended)) {
                            continue;
                        }
                    }
                    if (!IsMinimal(extended)) {
                        continue;
                    }
                    next_candidates[lhs].insert(extended);
                }
            }
        } else if (IsMinimal(lhs)) {
            AttributeList lhs_max_prefix = MaxPrefix(lhs);
            for (AttributeList const& rhs : candidate_sets_[lhs_max_prefix]) {
                if (AreDisjoint(lhs, rhs)) {
                    next_candidates[lhs].insert(rhs);
                }
            }
        }
        if (next_candidates[lhs].empty()) {
            next_candidates.erase(lhs);
        }
    }
    previous_candidate_sets_ = std::move(candidate_sets_);
    candidate_sets_ = std::move(next_candidates);
}

void Order::MergePrune() {
    if (lattice_->GetLevelNumber() < 3) {
        return;
    }
    for (auto const& [lhs, rhs_list] : candidate_sets_) {
        if (lhs.size() <= 1) {
            continue;
        }
        for (auto rhs_it = rhs_list.begin(); rhs_it != rhs_list.end();) {
            AttributeList lhs_max_prefix = MaxPrefix(lhs);
            if (InUnorderedMap(merge_invalidated_, lhs_max_prefix, *rhs_it)) {
                bool prunable = true;
                for (auto const& other_rhs : candidate_sets_[lhs_max_prefix]) {
                    if (MaxPrefix(other_rhs) == *rhs_it) {
                        prunable = false;
                        break;
                    }
                }
                if (prunable) {
                    rhs_it = candidate_sets_[lhs].erase(rhs_it);
                } else {
                    ++rhs_it;
                }
            } else {
                ++rhs_it;
            }
        }
    }
}

void Order::PrintValidOD() {
    LOG(DEBUG) << "***PREVIOUS CANDIDATE SETS***" << '\n';
    for (auto const& [lhs, rhs_list] : previous_candidate_sets_) {
        if (rhs_list.empty()) {
            for (auto const& attr : lhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << "-> empty";
            LOG(DEBUG) << '\n';
        }
        for (AttributeList const& rhs : rhs_list) {
            for (auto const& attr : lhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << "->";
            for (auto const& attr : rhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << '\n';
        }
    }
    LOG(DEBUG) << "***CANDIDATE SETS***" << '\n';
    for (auto const& [lhs, rhs_list] : candidate_sets_) {
        if (rhs_list.empty()) {
            for (auto const& attr : lhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << "-> empty";
            LOG(DEBUG) << '\n';
        }
        for (AttributeList const& rhs : rhs_list) {
            for (auto const& attr : lhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << "->";
            for (auto const& attr : rhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << '\n';
        }
    }
    LOG(DEBUG) << "***VALID ORDER DEPENDENCIES***" << '\n';
    unsigned int cnt = 0;
    for (auto const& [lhs, rhs_list] : valid_) {
        for (AttributeList const& rhs : rhs_list) {
            ++cnt;
            for (auto const& attr : lhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << "->";
            for (auto const& attr : rhs) {
                LOG(DEBUG) << attr + 1 << ",";
            }
            LOG(DEBUG) << '\n';
        }
    }
    LOG(DEBUG) << "OD amount: " << cnt;
    LOG(DEBUG) << '\n' << '\n';
}

unsigned long long Order::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();
    CreateSingleColumnSortedPartitions();
    lattice_ = std::make_unique<ListLattice>(candidate_sets_, single_attributes_);
    while (!lattice_->IsEmpty()) {
        ComputeDependencies(lattice_->GetLatticeLevel());
        lattice_->Prune(candidate_sets_);
        lattice_->GenerateNextLevel(candidate_sets_);
    }
    PrintValidOD();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    LOG(DEBUG) << "ms: " << elapsed_milliseconds.count() << '\n';
    return elapsed_milliseconds.count();
}

}  // namespace algos::order

#include "aid.h"

#include "config/tabular_data/input_table/option.h"

namespace algos {

Aid::Aid() : FDAlgorithm({kDefaultPhaseName}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName()});
}

void Aid::RegisterOptions() {
    RegisterOption(config::kTableOpt(&input_table_));
}

void Aid::LoadDataInternal() {
    number_of_attributes_ = input_table_->GetNumberOfColumns();
    if (number_of_attributes_ == 0) {
        throw std::runtime_error("Unable to work on an empty dataset.");
    }

    schema_ = RelationalSchema::CreateFrom(*input_table_);

    while (input_table_->HasNextRow()) {
        std::vector<std::string> const& next_line = input_table_->GetNextRow();
        if (next_line.empty()) {
            break;
        }

        tuples_.emplace_back(std::vector<size_t>(number_of_attributes_));
        for (size_t i = 0; i < number_of_attributes_; ++i) {
            tuples_.back()[i] = std::hash<std::string>{}(next_line[i]);
        }
    }
    number_of_tuples_ = tuples_.size();
    constant_columns_ = boost::dynamic_bitset<>(number_of_attributes_);
}

void Aid::ResetStateFd() {
    clusters_.assign(number_of_attributes_, std::unordered_map<size_t, Cluster>{});
    indices_in_clusters_.assign(number_of_attributes_, std::vector<size_t>(number_of_tuples_));
    constant_columns_.reset();
    prev_ratios_.assign(kWindowSize, 1.0);
    sum_ = double{kWindowSize};
}

unsigned long long Aid::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();

    BuildClusters();

    CreateNegativeCover();

    InvertNegativeCover();

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    return elapsed_milliseconds.count();
}

void Aid::BuildClusters() {
    if (number_of_tuples_ == 0) {
        return;
    }

    for (size_t attr_num = 0; attr_num < number_of_attributes_; ++attr_num) {
        std::unordered_map<size_t, Cluster>& column_values = clusters_[attr_num];
        bool is_constant = true;
        size_t first_value = tuples_[0][attr_num];
        for (size_t tuple_num = 0; tuple_num < number_of_tuples_; ++tuple_num) {
            size_t entry_value = tuples_[tuple_num][attr_num];
            if (entry_value != first_value) {
                is_constant = false;
            }

            Cluster& cluster = column_values[entry_value];
            cluster.push_back(tuple_num);
            indices_in_clusters_[attr_num][tuple_num] = cluster.size() - 1;
        }

        if (is_constant) {
            constant_columns_[attr_num] = true;
        }
    }
}

bool Aid::IsNegativeCoverGrowthSmall(size_t iteration_num, double curr_ratio) {
    iteration_num %= kWindowSize;
    sum_ -= prev_ratios_[iteration_num];
    sum_ += curr_ratio;
    prev_ratios_[iteration_num] = curr_ratio;
    double average = sum_ / (double)kWindowSize;
    if (average < kGrowthThreshold) {
        return true;
    }

    return false;
}

void Aid::CreateNegativeCover() {
    size_t prev_neg_cover_size = 0;
    for (size_t index = 1;; ++index) {
        for (size_t tuple_num = 0; tuple_num < number_of_tuples_; ++tuple_num) {
            HandleTuple(tuple_num, index);
        }

        size_t curr_neg_cover_size = neg_cover_.size();
        double curr_ratio;
        if (prev_neg_cover_size == 0) {
            curr_ratio = (curr_neg_cover_size == 0) ? 0.0 : 1.0;
        } else {
            curr_ratio = (double)curr_neg_cover_size / prev_neg_cover_size - 1;
        }

        // TODO: probably add other criteria of termination
        if (IsNegativeCoverGrowthSmall(index, curr_ratio)) {
            break;
        }

        prev_neg_cover_size = curr_neg_cover_size;
    }
}

void Aid::HandleTuple(size_t tuple_num, size_t iteration_num) {
    for (size_t attr_num = 0; attr_num < number_of_attributes_; ++attr_num) {
        size_t value = tuples_[tuple_num][attr_num];
        Cluster const& cluster = clusters_[attr_num].at(value);
        size_t index_in_cluster = indices_in_clusters_[attr_num][tuple_num];
        if (iteration_num <= index_in_cluster) {
            size_t another_index_in_cluster =
                    GenerateSecondClusterIndex(index_in_cluster, iteration_num);
            size_t another_tuple_num = cluster[another_index_in_cluster];
            auto tuples_agree_set = BuildAgreeSet(tuple_num, another_tuple_num);
            neg_cover_.insert(tuples_agree_set);
        }
    }
}

boost::dynamic_bitset<> Aid::BuildAgreeSet(size_t t1, size_t t2) {
    boost::dynamic_bitset<> equal_attr(number_of_attributes_);
    for (size_t attr_num = 0; attr_num < number_of_attributes_; ++attr_num) {
        if (tuples_[t1][attr_num] == tuples_[t2][attr_num]) {
            equal_attr.set(attr_num);
        }
    }

    return equal_attr;
}

void Aid::HandleConstantColumns(boost::dynamic_bitset<>& attributes) {
    boost::dynamic_bitset<> empty_set(number_of_attributes_);
    Vertical lhs = *schema_->empty_vertical_;
    for (size_t attr_num = constant_columns_.find_first();
         attr_num != boost::dynamic_bitset<>::npos;
         attr_num = constant_columns_.find_next(attr_num)) {
        attributes[attr_num] = false;
        Column rhs = schema_->GetColumn(attr_num);
        RegisterFd(lhs, rhs);
    }
}

void Aid::HandleInvalidFd(boost::dynamic_bitset<> const& neg_cover_el, SearchTree& pos_cover_tree,
                          size_t rhs) {
    std::vector<boost::dynamic_bitset<>> subsets;
    pos_cover_tree.ForEachSubset(neg_cover_el, [&subsets](boost::dynamic_bitset<> const& subset) {
        subsets.push_back(subset);
    });

    for (auto const& subset : subsets) {
        pos_cover_tree.Remove(subset);
    }

    for (auto& subset : subsets) {
        for (size_t add_attr = 0; add_attr < number_of_attributes_; ++add_attr) {
            if (add_attr == rhs || neg_cover_el[add_attr] || constant_columns_[add_attr]) {
                continue;
            }

            subset[add_attr] = true;
            if (!pos_cover_tree.ContainsAnySubsetOf(subset)) {
                pos_cover_tree.Add(subset);
            }

            subset[add_attr] = false;
        }
    }
}

void Aid::InvertNegativeCover() {
    boost::dynamic_bitset<> attributes(number_of_attributes_);
    attributes.set();
    HandleConstantColumns(attributes);

    std::vector<boost::dynamic_bitset<>> neg_cover_vector;
    neg_cover_vector.insert(neg_cover_vector.end(), neg_cover_.begin(), neg_cover_.end());
    auto comp_by_card = [](boost::dynamic_bitset<> const& lhs, boost::dynamic_bitset<> const& rhs) {
        return lhs.count() > rhs.count();
    };
    std::sort(neg_cover_vector.begin(), neg_cover_vector.end(), comp_by_card);
    auto attr_indices = GetAttributesSortedByFrequency(neg_cover_vector);

    std::vector<size_t> inv_attr_indices(number_of_attributes_);
    for (size_t i = 0; i < number_of_attributes_; ++i) {
        inv_attr_indices[attr_indices[i]] = i;
    }

    this->constant_columns_ = ChangeAttributesOrder(this->constant_columns_, inv_attr_indices);
    attributes = ChangeAttributesOrder(attributes, inv_attr_indices);
    for (auto& neg_cover_el : neg_cover_vector) {
        neg_cover_el = ChangeAttributesOrder(neg_cover_el, inv_attr_indices);
    }

    for (size_t rhs = 0; rhs < number_of_attributes_; ++rhs) {
        if (constant_columns_[rhs]) {
            continue;
        }

        attributes[rhs] = false;
        SearchTree pos_cover_tree(attributes);
        for (auto const& neg_cover_el : neg_cover_vector) {
            if (!neg_cover_el[rhs]) {
                HandleInvalidFd(neg_cover_el, pos_cover_tree, rhs);
            }
        }

        size_t real_rhs = attr_indices[rhs];
        std::vector<boost::dynamic_bitset<>> pos_cover_vector;
        pos_cover_tree.ForEach(
                [&pos_cover_vector, &attr_indices](boost::dynamic_bitset<> const& pos_cover_el) {
                    pos_cover_vector.push_back(ChangeAttributesOrder(pos_cover_el, attr_indices));
                });

        RegisterFDs(real_rhs, pos_cover_vector);
        attributes[rhs] = true;
    }
}

void Aid::RegisterFDs(size_t rhs_attribute,
                      std::vector<boost::dynamic_bitset<>> const& list_of_lhs_attributes) {
    Column rhs = schema_->GetColumn(rhs_attribute);
    for (auto const& lhs_attributes : list_of_lhs_attributes) {
        Vertical lhs = schema_->GetVertical(lhs_attributes);
        RegisterFd(lhs, rhs);
    }
}

size_t Aid::GenerateSecondClusterIndex(size_t index_in_cluster, size_t iteration_num) const {
    return (iteration_num * kPrime) % index_in_cluster;
}

boost::dynamic_bitset<> Aid::ChangeAttributesOrder(boost::dynamic_bitset<> const& initial_bitset,
                                                   std::vector<size_t> const& new_order) {
    size_t bitset_size = initial_bitset.size();
    boost::dynamic_bitset<> modified_bitset(bitset_size);
    for (size_t bit = 0; bit < bitset_size; ++bit) {
        if (initial_bitset[bit]) {
            modified_bitset.set(new_order[bit]);
        }
    }

    return modified_bitset;
}

std::vector<size_t> Aid::GetAttributesSortedByFrequency(
        std::vector<boost::dynamic_bitset<>> const& neg_cover_vector) const {
    std::vector<unsigned int> frequency(number_of_attributes_, 0);
    for (auto const& bitset : neg_cover_vector) {
        for (size_t attr_num = 0; attr_num < number_of_attributes_; ++attr_num) {
            frequency[attr_num] += bitset[attr_num];
        }
    }

    std::vector<size_t> attr_indices(number_of_attributes_);
    for (size_t attr_num = 0; attr_num < number_of_attributes_; ++attr_num) {
        attr_indices[attr_num] = attr_num;
    }

    std::sort(attr_indices.begin(), attr_indices.end(),
              [&frequency](size_t lhs, size_t rhs) { return frequency[lhs] > frequency[rhs]; });

    return attr_indices;
}

}  // namespace algos

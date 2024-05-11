#include "eulerfd.h"

#include <cstdlib>
#include <iostream>
#include <memory>
#include <limits>

namespace algos {

EulerFD::EulerFD() :
  FDAlgorithm({kDefaultPhaseName}), mlfq_(queues_number_) {
    last_ncover_ratios_.fill(1);
    last_pcover_ratios_.fill(1);
    RegisterOption(config::CustomRandomFlagOpt(&custom_random_opt_));

    RegisterOption(config::kTableOpt(&input_table_));
    MakeOptionsAvailable({config::kTableOpt.GetName()});

    max_lhs_ = std::numeric_limits<unsigned int>::max();
}

void EulerFD::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable({config::CustomRandomFlagOpt.GetName()});
}

void EulerFD::LoadDataInternal() {
    number_of_attributes_ = input_table_->GetNumberOfColumns();
    if (number_of_attributes_ == 0) {
        throw std::runtime_error("Unable to work on an empty dataset.");
    }

    schema_ = std::make_unique<RelationalSchema>(input_table_->GetRelationName());
    for (size_t i = 0; i < number_of_attributes_; ++i) {
        std::string const& column_name = input_table_->GetColumnName(static_cast<int>(i));
        schema_->AppendColumn(column_name);
    }

    // in each column mapping string values into integer values.
    // using only hash isnt good idea because colisions dont processing
    std::vector<std::unordered_map<std::string, size_t>> columns(number_of_attributes_);
    std::vector<size_t> current_id(number_of_attributes_, 0);

    while (input_table_->HasNextRow()) {
        const std::vector<std::string> line = input_table_->GetNextRow();
        if (line.empty()) {
            break;
        }

        tuples_.emplace_back(std::vector<size_t>(number_of_attributes_));
        for (size_t i = 0; i < number_of_attributes_; i++) {
            auto &values = columns[i];
            auto it = values.find(line[i]);
            if (it != values.end()) {
                tuples_.back()[i] = it->second;
            } else {
                size_t id = current_id[i]++;
                values[std::move(line[i])] = id;
                tuples_.back()[i] = id;
            }
        }
    }
    number_of_tuples_ = tuples_.size();
}

void EulerFD::ResetStateFd() {
    // data from sampling module
    clusters_.clear();
    constant_columns_.clear();
    is_first_sample_ = true;

    mlfq_.clear();
    effective_treshold_ = initial_effective_treshold_;

    invalids_.clear();
    new_invalids_.clear();
    fd_num_ = old_invalid_size_ = 0;

    // data from p/ncover consturction modules
    attribute_indexes_.clear();
    attribute_frequences_.clear();

    positive_cover_.clear();
    negative_cover_.clear();

    last_ncover_ratios_.fill(1);
    last_pcover_ratios_.fill(1);

    random_ = nullptr;
}

void EulerFD::SaveAnswer() {
    if (attribute_indexes_.size() == 0) {
        std::cout << "attribute_indexes_ size is 0\n";
        return;
    }

    std::vector<size_t> inv_indexes(number_of_attributes_);
    for (size_t i = 0; i < number_of_attributes_; ++i) {
        inv_indexes[attribute_indexes_[i]] = i;
    }

    for (size_t rhs_attr = 0; rhs_attr < number_of_attributes_; rhs_attr++) {
        // then tree filling we use inverse indexes, so to get correct tree we inverse again
        size_t inv_rhs_attr = inv_indexes[rhs_attr];
        auto &tree = positive_cover_[inv_rhs_attr];
        auto rhs = *schema_->GetColumn(rhs_attr);
        tree.ForEach(
            [&](const Bitset &lhs_attr) {
                auto lhs = schema_->GetVertical(ChangeAttributesOrder(lhs_attr, attribute_indexes_));
                RegisterFd(lhs, rhs);
            });
    }
}

void EulerFD::InitCovers() {
    positive_cover_.clear();
    positive_cover_.reserve(number_of_attributes_);
    for (size_t i = 0; i < number_of_attributes_; i++) {
        positive_cover_.push_back(SearchTreeEulerFD(number_of_attributes_));
        positive_cover_.back().Add(Bitset(number_of_attributes_));
    }

    negative_cover_.clear();
    negative_cover_.reserve(number_of_attributes_);
    for (size_t i = 0; i < number_of_attributes_; i++) {
        negative_cover_.push_back(SearchTreeEulerFD(number_of_attributes_));
    }

    attribute_frequences_.resize(number_of_attributes_);
    std::fill(attribute_frequences_.begin(), attribute_frequences_.end(), 0);

    attribute_indexes_.clear();
}

void EulerFD::BuildPartition() {
    if (number_of_tuples_ == 0) {
        return;
    }
    constant_columns_ = Bitset(number_of_attributes_);

    // build partition: same values in each column unite into clusters,
    // then stripping it: clusters with 1 element cant get us new invalid fd
    for (size_t attr_num = 0; attr_num < number_of_attributes_; attr_num++) {
        std::unordered_map<size_t, std::vector<size_t>> values;
        for (size_t tuple_num = 0; tuple_num < number_of_tuples_; tuple_num++) {
            auto value = tuples_[tuple_num][attr_num];
            auto &similar_values = values[value];
            similar_values.push_back(tuple_num);
        }
        if (values.size() == 1) {
            // its mean in this column only one value
            constant_columns_[attr_num] = true;
            continue;
        } else if (values.size() == number_of_tuples_) {
            // its mean all clusters' sizes are 1, so we can pass cycle below
            continue;
        } else {
            for (auto &&[_, cluster] : values) {
                if (cluster.size() > 1) {
                    clusters_.emplace_back(std::move(cluster), rand_function_);
                }
            }
        }
    }
}

EulerFD::Bitset EulerFD::BuildAgreeSet(size_t t1, size_t t2) const {
    Bitset equal_attr(number_of_attributes_);
    for (size_t attr_num = 0; attr_num < number_of_attributes_; ++attr_num) {
        if (tuples_[t1][attr_num] == tuples_[t2][attr_num]) {
            equal_attr.set(attr_num);
        }
    }
    return equal_attr;
}

double EulerFD::SamlingInCluster(Cluster *cluster) {
    return cluster->Sample(
        [this](size_t t1, size_t t2) -> size_t {
            Bitset agree_set = BuildAgreeSet(t1, t2);
            auto &&[_, result] = invalids_.insert(agree_set);

            // check is discovered fd new
            if (result) {
                new_invalids_.insert(agree_set);
                return agree_set.size() - agree_set.count();
            } else {
                return 0;
            }
        });
}

void EulerFD::Sampling() {
    if (is_first_sample_) {
        // in first sampling mlfq is empty, we fiil it.
        // ww put all clusters in mlfq, even if effective was 0
        is_first_sample_ = false;

        for (size_t i = 0; i < clusters_.size(); i++) { // TODO: maybe ranged for
            double eff = SamlingInCluster(&clusters_[i]);
            mlfq_.Add(&clusters_[i], eff, true);
        }

        if (mlfq_.GetLastQueueSize() > 0) {
            effective_treshold_ = std::min(effective_treshold_,
                mlfq_.MaxEffectInLastQueue() / 2);
        }

        return;
    }

    // sampling in first queues of mlfq;
    new_invalids_.clear();
    while (mlfq_.GetEffectiveSize() > 0) {
       Cluster *cluster = mlfq_.Get();
       SamlingInCluster(cluster);
       mlfq_.Add(cluster, cluster->GetAverage());
    }

    if (mlfq_.GetLastQueueSize() > 0) {
        effective_treshold_ = std::min(effective_treshold_ / 4,
            mlfq_.MaxEffectInLastQueue() / 2);
    }

    // sampling in last queues (it is priority queue) of mlfq;
    while (mlfq_.GetLastQueueSize() > 0 &&
      mlfq_.MaxEffectInLastQueue() >= effective_treshold_) {
        Cluster *cluster = mlfq_.Get();
        SamlingInCluster(cluster);
        mlfq_.AddAtLast(cluster);
    }
}

bool EulerFD::IsNCoverGrowthSmall() const {
    double sum = std::accumulate(last_ncover_ratios_.begin(), last_ncover_ratios_.end(), 0.0);
    double grow_rate = sum / window_;
    return grow_rate < neg_cover_growth_treshold_;
}

bool EulerFD::IsPCoverGrowthSmall() const {
    double sum = std::accumulate(last_pcover_ratios_.begin(), last_pcover_ratios_.end(), 0.0);
    double grow_rate = sum / window_;
    return grow_rate < pos_cover_growth_treshold_;
}

std::vector<size_t> EulerFD::GetAttributesSortedByFrequency(
  const std::vector<Bitset> &neg_cover_vector) {

    for (auto const& bitset : neg_cover_vector) {
        for (size_t attr_num = 0; attr_num < number_of_attributes_; ++attr_num) {
            attribute_frequences_[attr_num] += bitset[attr_num];
        }
    }

    std::vector<size_t> attr_indices(number_of_attributes_);
    std::iota(attr_indices.begin(), attr_indices.end(), 0);

    // after sorting in indexes[0] is number of column with the largest number of 1
    std::sort(attr_indices.begin(), attr_indices.end(),
        [&](size_t lhs, size_t rhs) {
            return attribute_frequences_[lhs] < attribute_frequences_[rhs];
        });

    return attr_indices;
}

EulerFD::Bitset EulerFD::ChangeAttributesOrder(
  const Bitset &initial_bitset, const std::vector<size_t> &new_order) {
    size_t bitset_size = initial_bitset.size();
    Bitset modified_bitset(bitset_size);
    for (size_t bit = 0; bit < bitset_size; ++bit) {
        if (initial_bitset[bit]) {
            modified_bitset.set(new_order[bit]);
        }
    }
    return modified_bitset;
}

void EulerFD::AddInvalidAtTree(SearchTreeEulerFD &tree, const Bitset &invalid) {
    if (tree.ContainsAnySupersetOf(invalid)) {
        return;
    }
    RemoveGeneralizations(tree, invalid);
    tree.Add(invalid);
}

std::unordered_set<EulerFD::Bitset> EulerFD::RemoveGeneralizations(
  SearchTreeEulerFD &tree, const Bitset &invalid) {
    // Generalizations is subsets
    std::unordered_set<Bitset> remove;
    tree.ForEachSubset(invalid,
        [&remove](const Bitset &sub) {
             remove.insert(sub);
        });
    for (auto &sub : remove) {
        tree.Remove(sub);
    }
    return remove;
}

std::vector<EulerFD::Bitset> EulerFD::CreateNegativeCover(size_t rhs,
  const std::vector<Bitset> &neg_cover_vector) {
    auto &tree = negative_cover_[rhs];
    for (auto &invalid : neg_cover_vector) {
        if (!invalid[rhs]) {
            AddInvalidAtTree(tree, invalid);
        }
    }

    std::vector<Bitset> neg_cover_rhs;
    tree.ForEach(
        [&neg_cover_rhs](const Bitset &invalid) {
            neg_cover_rhs.push_back(invalid);
        });
    return neg_cover_rhs;
}

size_t EulerFD::Invert(size_t rhs, const std::vector<Bitset> &neg) {
    auto &tree = positive_cover_[rhs];
    for (auto &invalid : neg) {
        auto removeds = RemoveGeneralizations(tree, invalid);
        for (auto &removed : removeds) {
            for (size_t i = 0; i < number_of_attributes_; i++) {
                if (i == rhs || invalid[i] || constant_columns_[attribute_indexes_[i]]) {
                    continue;
                }
                Bitset copy = removed;
                copy.set(i);
                if (!tree.ContainsAnySubsetOf(copy)) {
                    tree.Add(copy);
                }
            }
        }
    }
    return tree.GetCardinality();
}

size_t EulerFD::GenerateResults() {
    // check is new non fd discovered
    if (old_invalid_size_ == invalids_.size()) {
        return fd_num_;
    }
    old_invalid_size_ = invalids_.size();

    std::vector<Bitset> neg_cover_vector(new_invalids_.begin(), new_invalids_.end());
    auto new_indexes = GetAttributesSortedByFrequency(neg_cover_vector);

    // check is rebuild of covers neccesary
    if (attribute_indexes_ != new_indexes) {
        if (old_invalid_size_ != 0) {
            InitCovers();
            neg_cover_vector.assign(invalids_.begin(), invalids_.end());
            new_indexes = GetAttributesSortedByFrequency(neg_cover_vector);
        }
        attribute_indexes_ = new_indexes;
    }

    std::vector<size_t> inv_indexes(number_of_attributes_);
    for (size_t i = 0; i < number_of_attributes_; ++i) {
        inv_indexes[attribute_indexes_[i]] = i;
    }

    // change attribute order, in first index the largset number of 1
    for (auto& neg_cover_el : neg_cover_vector) {
        neg_cover_el = ChangeAttributesOrder(neg_cover_el, inv_indexes);
    }

    // sorting all non fd by number of 1
    std::sort(neg_cover_vector.begin(), neg_cover_vector.end(),
        [](const Bitset &left, const Bitset &right) {
            return left.count() > right.count();
        });

    // creating ncover and pcover trees for each rhs
    size_t fd_num = 0;
    for (size_t rhs = 0; rhs < number_of_attributes_; rhs++) {
        if (constant_columns_[rhs]) {
            continue;
        }

        size_t real_rhs = inv_indexes[rhs];
        auto neg = CreateNegativeCover(real_rhs, neg_cover_vector);
        // sorting all non fd by number of 1
        std::sort(neg.begin(), neg.end(),
            [](const Bitset &left, const Bitset &right) {
                return left.count() > right.count();
            });
        fd_num += Invert(real_rhs, neg);
    }
    return fd_num;
}

unsigned long long EulerFD::ExecuteInternal() {
    if (number_of_attributes_ == 1) {
        return 0;
    }

    // choose random strategy (it is neccesary for stable unit tests)
    if (custom_random_opt_.first) {
        random_ = std::make_unique<CustomRandom>(custom_random_opt_.second);
        rand_function_ = [&](){ return random_->NextInt(random_upper_bound_); };
    } else {
        srand(time(NULL));
        rand_function_ = std::rand;
    }

    auto start_time = std::chrono::system_clock::now();

    BuildPartition();
    if (clusters_.size() == 0) {
        // in small datasets sometimes after clusters stripping there are no clusters for sampling
        std::cout << "number of clusters is 0*\n";
        return 0;
    }

    InitCovers();
    size_t iteration_number = 0;
    while (true) {
        size_t ncover_size = invalids_.size();
        Sampling();

        last_ncover_ratios_[iteration_number % window_] =
            invalids_.size() == 0 ? 0 :
            (double)(invalids_.size() - ncover_size) / invalids_.size();

        // check criterion for enter second EulerFD cycle
        if (IsNCoverGrowthSmall()) {
            size_t pcover_size = fd_num_;
            fd_num_ = GenerateResults();
            last_pcover_ratios_[iteration_number % window_] =
                fd_num_ == 0 ? 0 :
                (double)(fd_num_ - pcover_size) / fd_num_;
            if (IsPCoverGrowthSmall()) {
                break;
            }
        }

        iteration_number++;
    }

    // Convert answer from pcover trees to list of fd
    SaveAnswer();

    auto elapsed_milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}
} // end of 'algos' namespace

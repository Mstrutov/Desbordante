//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <map>
#include <memory>
#include <utility>

#include <boost/dynamic_bitset.hpp>
#include <easylogging++.h>

#include "ColumnLayoutRelationData.h"
#include "PositionListIndex.h"
#include "Vertical.h"

namespace util {

const int PositionListIndex::singleton_value_id_ = 0;
unsigned long long PositionListIndex::micros_ = 0;
int PositionListIndex::intersection_count_ = 0;

PositionListIndex::PositionListIndex(std::deque<std::vector<int>> index,
                                     std::vector<int> null_cluster, unsigned int size,
                                     double entropy, unsigned long long nep,
                                     unsigned int relation_size,
                                     unsigned int original_relation_size, double inverted_entropy,
                                     double gini_impurity)
    : index_(std::move(index)),
      null_cluster_(std::move(null_cluster)),
      size_(size),
      entropy_(entropy),
      inverted_entropy_(inverted_entropy),
      gini_impurity_(gini_impurity),
      nep_(nep),
      relation_size_(relation_size),
      original_relation_size_(original_relation_size),
      probing_table_cache_() {}

std::unique_ptr<PositionListIndex> PositionListIndex::CreateFor(std::vector<int>& data,
                                                                bool is_null_eq_null) {
    std::unordered_map<int, std::vector<int>> index;
    for (unsigned long position = 0; position < data.size(); ++position) {
        int value_id = data[position];
        index[value_id].push_back(position);
    }

    std::vector<int> null_cluster;
    if (index.count(ColumnLayoutRelationData::kNullValueId) != 0) {
        null_cluster = index[ColumnLayoutRelationData::kNullValueId];
    }
    if (!is_null_eq_null) {
        index.erase(ColumnLayoutRelationData::kNullValueId); // move?
    }

    double key_gap = 0.0;
    double inv_ent = 0;
    double gini_gap = 0;
    unsigned long long nep = 0;
    unsigned int size = 0;
    std::deque<std::vector<int>> clusters;

    for (auto& iter : index) {
        if (iter.second.size() == 1) {
            gini_gap += std::pow(1 / static_cast<double>(data.size()), 2);
            continue;
        }
        key_gap += iter.second.size() * log(iter.second.size());
        nep += CalculateNep(iter.second.size());
        size += iter.second.size();
        inv_ent += -(1 - iter.second.size() / static_cast<double>(data.size())) *
                   std::log(1 - (iter.second.size() / static_cast<double>(data.size())));
        gini_gap += std::pow(iter.second.size() / static_cast<double>(data.size()), 2);

        clusters.emplace_back(std::move(iter.second));
    }
    double entropy = log(data.size()) - key_gap / data.size();

    double gini_impurity = 1 - gini_gap;
    if (gini_impurity == 0) {
        inv_ent = 0;
    }

    SortClusters(clusters);
    return std::make_unique<PositionListIndex>(std::move(clusters), std::move(null_cluster), size,
                                               entropy, nep, data.size(), data.size(), inv_ent,
                                               gini_impurity);
}

//unsigned long long PositionListIndex::CalculateNep(unsigned int numElements) {
//
//}

void PositionListIndex::SortClusters(std::deque<std::vector<int>>& clusters) {
    sort(clusters.begin(), clusters.end(),
         [](std::vector<int>& a, std::vector<int>& b) { return a[0] < b[0]; });
}

std::shared_ptr<const std::vector<int>> PositionListIndex::CalculateAndGetProbingTable() const {
    if (probing_table_cache_ != nullptr) return probing_table_cache_;

    std::vector<int> probing_table = std::vector<int>(original_relation_size_);
    int next_cluster_id = singleton_value_id_ + 1;
    for (auto& cluster : index_) {
        int value_id = next_cluster_id++;
        assert(value_id != singleton_value_id_);
        for (int position : cluster) {
            probing_table[position] = value_id;
        }
    }

    return std::make_shared<std::vector<int>>(probing_table);
}



// интересное место: true --> надо передать поле без копирования, false --> надо сконструировать и выдать наружу
// кажется, самым лёгким способом будет навернуть shared_ptr
/*std::shared_ptr<const std::vector<int>> PositionListIndex::getProbingTable(bool isCaching) {
    auto probingTable = GetProbingTable();
    if (isCaching) {
        probing_table_cache_ = probingTable;
        return probing_table_cache_;
    }
    return probingTable;
}*/

//std::deque<std::vector<int>> const & PositionListIndex::getIndex() const {
//    return index;
//}

std::unique_ptr<PositionListIndex> PositionListIndex::Intersect(PositionListIndex const* that) const {
    assert(this->relation_size_ == that->relation_size_);
    return this->size_ > that->size_ ?
           that->Probe(this->CalculateAndGetProbingTable()) :
           this->Probe(that->CalculateAndGetProbingTable());
}

//TODO: null_cluster_ некорректен
std::unique_ptr<PositionListIndex> PositionListIndex::Probe(std::shared_ptr<const std::vector<int>> probing_table) const {
    assert(this->relation_size_ == probing_table->size());
    std::deque<std::vector<int>> new_index;
    unsigned int new_size = 0;
    double new_key_gap = 0.0;
    unsigned long long new_nep = 0;
    std::vector<int> null_cluster;

    std::unordered_map<int, std::vector<int>> partial_index;


    for (auto& positions : index_) {
        for (int position : positions) {
            if (probing_table == nullptr) LOG(DEBUG) << "NULLPTR";
            if (position < 0 || static_cast<size_t>(position) >= probing_table->size()) {
                LOG(DEBUG) << "position: " + std::to_string(position) +
                                 ", size: " + std::to_string(probing_table->size());
                for (size_t i = 0; i < positions.size(); ++i) {
                    LOG(DEBUG) << "Position " + std::to_string(positions[i]);
                }
            }
            int probing_table_value_id = (*probing_table)[position];
            if (probing_table_value_id == singleton_value_id_)
                continue;
            intersection_count_++;
            partial_index[probing_table_value_id].push_back(position);
        }

        for (auto& iter : partial_index) {
            auto& cluster = iter.second;
            if (cluster.size() <= 1) continue;

            new_size += cluster.size();
            new_key_gap += cluster.size() * log(cluster.size());
            new_nep += CalculateNep(cluster.size());

            new_index.push_back(std::move(cluster));
        }
        partial_index.clear();
    }

    double new_entropy = log(relation_size_) - new_key_gap / relation_size_;
    SortClusters(new_index);

    return std::make_unique<PositionListIndex>(std::move(new_index), std::move(null_cluster),
                                               new_size, new_entropy, new_nep, relation_size_,
                                               relation_size_);
}

//TODO: null_cluster_ не поддерживается
std::unique_ptr<PositionListIndex> PositionListIndex::ProbeAll(
    Vertical const& probing_columns, ColumnLayoutRelationData& relation_data) {
    assert(this->relation_size_ == relation_data.GetNumRows());
    std::deque<std::vector<int>> new_index;
    unsigned int new_size = 0;
    double new_key_gap = 0.0;
    unsigned long long new_nep = 0;

    std::map<std::vector<int>, std::vector<int>> partial_index;
    std::vector<int> null_cluster;
    std::vector<int> probe;

    for (auto& cluster : this->index_) {
        for (int position : cluster) {
            if (!TakeProbe(position, relation_data, probing_columns, probe)) {
                probe.clear();
                continue;
            }

            partial_index[probe].push_back(position);
            probe.clear();
        }

        for (auto& iter : partial_index) {
            auto& new_cluster = iter.second;
            if (new_cluster.size() == 1) continue;

            new_size += new_cluster.size();
            new_key_gap += new_cluster.size() * log(new_cluster.size());
            new_nep += CalculateNep(new_cluster.size());

            new_index.emplace_back(std::move(new_cluster));
        }
        partial_index.clear();
    }

    double new_entropy = log(this->relation_size_) - new_key_gap / this->relation_size_;

    SortClusters(new_index);

    return std::make_unique<PositionListIndex>(std::move(new_index), std::move(null_cluster),
                                               new_size, new_entropy, new_nep, this->relation_size_,
                                               this->relation_size_);
}

bool PositionListIndex::TakeProbe(int position, ColumnLayoutRelationData& relation_data,
                                  Vertical const& probing_columns, std::vector<int>& probe) {
    boost::dynamic_bitset<> probing_indices = probing_columns.GetColumnIndices();
    for (unsigned long index = probing_indices.find_first(); index < probing_indices.size();
         index = probing_indices.find_next(index)) {
        int value = relation_data.GetColumnData(index).GetProbingTableValue(position);
        if (value == PositionListIndex::singleton_value_id_) return false;
        probe.push_back(value);
    }
    return true;
}

std::string PositionListIndex::ToString() const {
    std::string res = "[";
    for (auto& cluster : index_) {
        res.push_back('[');
        for (int v : cluster) {
            res.append(std::to_string(v) + ",");
        }
        if (res.find(',') != std::string::npos)
            res.erase(res.find_last_of(','));
        res.push_back(']');
        res.push_back(',');
    }
    if (res.find(',') != std::string::npos)
        res.erase(res.find_last_of(','));
    res.push_back(']');
    return res;
}

} // namespace util


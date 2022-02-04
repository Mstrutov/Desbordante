#pragma once

#include <functional>
#include <unordered_map>
#include <random>

#include "AgreeSetSample.h"

#include "easylogging++.h"

namespace util {

template <typename T>
std::unique_ptr<T> AgreeSetSample::CreateFor(ColumnLayoutRelationData* relation_data, int sample_size) {
    static_assert(std::is_base_of<AgreeSetSample, T>::value);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> random(0, relation_data->GetNumRows());

    std::unordered_map<boost::dynamic_bitset<>, int> agree_set_counters;
  sample_size = std::min((unsigned long long)sample_size, relation_data->GetNumTuplePairs());

    for (long i = 0; i < sample_size; i++){
        int tuple_index_1 = random(gen);
        int tuple_index_2 = random(gen);
        if (tuple_index_1 == tuple_index_2){
            i--;
            continue;
        }

        boost::dynamic_bitset<> agree_set(relation_data->GetNumColumns());
        for (auto & column_data : relation_data->GetColumnData()){
            int value1 = column_data.GetProbingTableValue(tuple_index_1);
            if (value1 != PositionListIndex::singleton_value_id_ && value1 == column_data.GetProbingTableValue(tuple_index_2)){
                agree_set[column_data.GetColumn()->GetIndex()] = true;
            }
        }

        agree_set_counters[agree_set]++;
    }

    auto instance = std::make_unique<T>(relation_data, relation_data->GetSchema()->empty_vertical_,
                                        (int) sample_size, relation_data->GetNumTuplePairs(), agree_set_counters);
    return instance;
}

template<typename T>
std::unique_ptr<T> AgreeSetSample::CreateFocusedFor(ColumnLayoutRelationData const* relation,
                                                    Vertical const& restriction_vertical,
                                                    PositionListIndex const* restriction_pli, unsigned int sample_size,
                                                    CustomRandom& random) {
    static_assert(std::is_base_of<AgreeSetSample, T>::value);
    //std::random_device rd;
    //std::mt19937 gen(rd());
    //std::uniform_real_distribution<> random_double;

    boost::dynamic_bitset<> free_column_indices(relation->GetNumColumns());
    free_column_indices.set();
    free_column_indices &= ~restriction_vertical.GetColumnIndices();
    std::vector<std::reference_wrapper<const ColumnData>> relevant_column_data;
    for (size_t column_index = free_column_indices.find_first();
         column_index != boost::dynamic_bitset<>::npos;
         column_index = free_column_indices.find_next(column_index)){
        relevant_column_data.emplace_back(relation->GetColumnData(column_index));
    }
    boost::dynamic_bitset<> agree_set_prototype(restriction_vertical.GetColumnIndices());
    std::unordered_map<boost::dynamic_bitset<>, int> agree_set_counters;

    unsigned long long restriction_nep = restriction_pli->GetNepAsLong();
  sample_size = std::min(static_cast<unsigned long long>(sample_size), restriction_nep);
    if (sample_size >= restriction_nep){
        for (auto & cluster : restriction_pli->GetIndex()){
            for (unsigned int i = 0; i < cluster.size(); i++){
                int tuple_index_1 = cluster[i];
                for (unsigned int j = i + 1; j < cluster.size(); j++){
                    int tuple_index_2 = cluster[j];

                    boost::dynamic_bitset<> agree_set(agree_set_prototype);
                    for (auto & column_data : relevant_column_data){
                        int value1 = column_data.get().GetProbingTableValue(tuple_index_1);
                        if (value1 != PositionListIndex::singleton_value_id_
                             && value1 == column_data.get().GetProbingTableValue(tuple_index_2)){
                            agree_set.set(column_data.get().GetColumn()->GetIndex());
                        }
                    }
                    auto location = agree_set_counters.find(agree_set);
                    if (location == agree_set_counters.end()) {
                        agree_set_counters.emplace_hint(location, agree_set, 1);
                    } else {
                        location->second += 1;
                    }
                }
            }
        }
    } else {
        std::vector<unsigned long long> cluster_sizes(restriction_pli->GetNumNonSingletonCluster() - 1);
        for (unsigned int i = 0; i < cluster_sizes.size(); i++){
            unsigned long long cluster_size = restriction_pli->GetIndex()[i].size();
            unsigned long long num_tuple_pairs = cluster_size * (cluster_size - 1) / 2;
            if (i > 0){
                cluster_sizes[i] = num_tuple_pairs + cluster_sizes[i - 1];
            } else {
                cluster_sizes[i] = num_tuple_pairs;
            }
        }

        for (unsigned int i = 0; i < sample_size; i++){
            auto cluster_index_iter = std::lower_bound(cluster_sizes.begin(), cluster_sizes.end(),
                                                       random.NextULL() % restriction_nep);
            unsigned int cluster_index = std::distance(cluster_sizes.begin(), cluster_index_iter);
            /*if (cluster_index >= cluster_sizes.size()) {
                cluster_index = cluster_sizes.size() - 1;
            }*/
            auto& cluster = restriction_pli->GetIndex()[cluster_index];

            int tuple_index_1 = random.NextInt(cluster.size());
            int tuple_index_2 = random.NextInt(cluster.size());
            while (tuple_index_1 == tuple_index_2) {
                tuple_index_2 = random.NextInt(cluster.size());
            }
            tuple_index_1 = cluster[tuple_index_1];
            tuple_index_2 = cluster[tuple_index_2];


            boost::dynamic_bitset<> agree_set(agree_set_prototype);
            for (auto & column_data : relevant_column_data){
                int value1 = column_data.get().GetProbingTableValue(tuple_index_1);
                if (value1 != PositionListIndex::singleton_value_id_
                     && value1 == column_data.get().GetProbingTableValue(tuple_index_2)){
                    agree_set.set(column_data.get().GetColumn()->GetIndex());
                }
            }

            auto location = agree_set_counters.find(agree_set);
            if (location == agree_set_counters.end()) {
                agree_set_counters.emplace_hint(location, agree_set, 1);
            } else {
                location->second += 1;
            }
        }
    }
    //std::cout << "-----------------\n";
    /*string agreeSetCountersStr = "{";
    for (auto& [key, value] : agree_set_counters) {
        agreeSetCountersStr += '\"';
        for (unsigned int columnIndex = key.find_first(); columnIndex < key.size(); columnIndex = key.find_next(columnIndex)){
            agreeSetCountersStr += std::to_string(columnIndex) + ' ';
        }
        agreeSetCountersStr += '\"';
        agreeSetCountersStr += " : "+ std::to_string(value) + ',';
    }
    agreeSetCountersStr.erase(agreeSetCountersStr.end()-1);
    agreeSetCountersStr += '}';

    LOG(DEBUG) << boost::format {"Created sample focused on %1%: %2%"} % restriction_vertical->ToString() % agreeSetCountersStr;
    */
    return std::make_unique<T>(relation, restriction_vertical, sample_size, restriction_nep, std::move(agree_set_counters));
}

} // namespace util


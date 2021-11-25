#pragma once

#include <functional>
#include <unordered_map>
#include <random>

#include "AgreeSetSample.h"

#include "easylogging++.h"

namespace util {

template <typename T>
std::unique_ptr<T> AgreeSetSample::createFor(ColumnLayoutRelationData* relationData, int sampleSize) {
    static_assert(std::is_base_of<AgreeSetSample, T>::value);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> random(0, relationData->getNumRows());

    std::unordered_map<boost::dynamic_bitset<>, int> agreeSetCounters;
    sampleSize = std::min((unsigned long long)sampleSize, relationData->getNumTuplePairs());

    for (long i = 0; i < sampleSize; i++){
        int tupleIndex1 = random(gen);
        int tupleIndex2 = random(gen);
        if (tupleIndex1 == tupleIndex2){
            i--;
            continue;
        }

        boost::dynamic_bitset<> agreeSet(relationData->getNumColumns());
        for (auto & columnData : relationData->getColumnData()){
            int value1 = columnData.getProbingTableValue(tupleIndex1);
            if (value1 != PositionListIndex::singletonValueId && value1 == columnData.getProbingTableValue(tupleIndex2)){
                agreeSet[columnData.getColumn()->getIndex()] = true;
            }
        }

        agreeSetCounters[agreeSet]++;
    }

    auto instance = std::make_unique<T>(relationData, relationData->getSchema()->emptyVertical,
                                        (int) sampleSize, relationData->getNumTuplePairs(), agreeSetCounters);
    return instance;
}

template<typename T>
std::unique_ptr<T> AgreeSetSample::createFocusedFor(ColumnLayoutRelationData const* relation,
                                               Vertical const& restrictionVertical,
                                               PositionListIndex const* restrictionPli, unsigned int sampleSize,
                                               CustomRandom& random) {
    static_assert(std::is_base_of<AgreeSetSample, T>::value);
    //std::random_device rd;
    //std::mt19937 gen(rd());
    //std::uniform_real_distribution<> random_double;

    boost::dynamic_bitset<> freeColumnIndices(relation->getNumColumns());
    freeColumnIndices.set();
    freeColumnIndices &= ~restrictionVertical.getColumnIndices();
    std::vector<std::reference_wrapper<const ColumnData>> relevantColumnData;
    for (size_t columnIndex = freeColumnIndices.find_first();
         columnIndex != boost::dynamic_bitset<>::npos;
         columnIndex = freeColumnIndices.find_next(columnIndex)){
        relevantColumnData.emplace_back(relation->getColumnData(columnIndex));
    }
    boost::dynamic_bitset<> agreeSetPrototype(restrictionVertical.getColumnIndices());
    std::unordered_map<boost::dynamic_bitset<>, int> agreeSetCounters;

    unsigned long long restrictionNep = restrictionPli->getNepAsLong();
    sampleSize = std::min(static_cast<unsigned long long>(sampleSize), restrictionNep);
    if (sampleSize >= restrictionNep){
        for (auto & cluster : restrictionPli->getIndex()){
            for (unsigned int i = 0; i < cluster.size(); i++){
                int tupleIndex1 = cluster[i];
                for (unsigned int j = i + 1; j < cluster.size(); j++){
                    int tupleIndex2 = cluster[j];

                    boost::dynamic_bitset<> agreeSet(agreeSetPrototype);
                    for (auto & columnData : relevantColumnData){
                        int value1 = columnData.get().getProbingTableValue(tupleIndex1);
                        if (value1 != PositionListIndex::singletonValueId
                             && value1 == columnData.get().getProbingTableValue(tupleIndex2)){
                            agreeSet.set(columnData.get().getColumn()->getIndex());
                        }
                    }
                    auto location = agreeSetCounters.find(agreeSet);
                    if (location == agreeSetCounters.end()) {
                        agreeSetCounters.emplace_hint(location, agreeSet, 1);
                    } else {
                        location->second += 1;
                    }
                }
            }
        }
    } else {
        std::vector<unsigned long long> clusterSizes(restrictionPli->getNumNonSingletonCluster() - 1);
        for (unsigned int i = 0; i < clusterSizes.size(); i++){
            unsigned long long clusterSize = restrictionPli->getIndex()[i].size();
            unsigned long long numTuplePairs = clusterSize * (clusterSize - 1) / 2;
            if (i > 0){
                clusterSizes[i] = numTuplePairs + clusterSizes[i - 1];
            } else {
                clusterSizes[i] = numTuplePairs;
            }
        }

        for (unsigned int i = 0; i < sampleSize; i++){
            auto clusterIndexIter = std::lower_bound(clusterSizes.begin(), clusterSizes.end(),
                                                     random.nextULL() % restrictionNep);
            unsigned int clusterIndex = std::distance(clusterSizes.begin(), clusterIndexIter);
            /*if (clusterIndex >= clusterSizes.size()) {
                clusterIndex = clusterSizes.size() - 1;
            }*/
            auto& cluster = restrictionPli->getIndex()[clusterIndex];

            int tupleIndex1 = random.nextInt(cluster.size());
            int tupleIndex2 = random.nextInt(cluster.size());
            while (tupleIndex1 == tupleIndex2) {
                tupleIndex2 = random.nextInt(cluster.size());
            }
            tupleIndex1 = cluster[tupleIndex1];
            tupleIndex2 = cluster[tupleIndex2];


            boost::dynamic_bitset<> agreeSet(agreeSetPrototype);
            for (auto & columnData : relevantColumnData){
                int value1 = columnData.get().getProbingTableValue(tupleIndex1);
                if (value1 != PositionListIndex::singletonValueId
                     && value1 == columnData.get().getProbingTableValue(tupleIndex2)){
                    agreeSet.set(columnData.get().getColumn()->getIndex());
                }
            }

            auto location = agreeSetCounters.find(agreeSet);
            if (location == agreeSetCounters.end()) {
                agreeSetCounters.emplace_hint(location, agreeSet, 1);
            } else {
                location->second += 1;
            }
        }
    }
    //std::cout << "-----------------\n";
    /*string agreeSetCountersStr = "{";
    for (auto& [key, value] : agreeSetCounters) {
        agreeSetCountersStr += '\"';
        for (unsigned int columnIndex = key.find_first(); columnIndex < key.size(); columnIndex = key.find_next(columnIndex)){
            agreeSetCountersStr += std::to_string(columnIndex) + ' ';
        }
        agreeSetCountersStr += '\"';
        agreeSetCountersStr += " : "+ std::to_string(value) + ',';
    }
    agreeSetCountersStr.erase(agreeSetCountersStr.end()-1);
    agreeSetCountersStr += '}';

    LOG(DEBUG) << boost::format {"Created sample focused on %1%: %2%"} % restrictionVertical->toString() % agreeSetCountersStr;
    */
    return std::make_unique<T>(relation, restrictionVertical, sampleSize, restrictionNep, std::move(agreeSetCounters));
}

} // namespace util


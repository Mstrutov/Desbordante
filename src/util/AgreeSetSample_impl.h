#pragma once

#include <iterator>
#include <random>

#include "AgreeSetSample.h"
//#include <utility>
//TODO: рандом проверь
//template<typename T, typename enable_if<is_base_of<AgreeSetSample, T>::value>::type*>
template <typename T>
shared_ptr<T> AgreeSetSample::createFor(shared_ptr<ColumnLayoutRelationData> relationData, int sampleSize) {
    static_assert(std::is_base_of<AgreeSetSample, T>::value);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> random(0, relationData->getNumRows());

    map<dynamic_bitset<>, int> agreeSetCounters;
    sampleSize = std::min((unsigned long long)sampleSize, relationData->getNumTuplePairs());

    for (long i = 0; i < sampleSize; i++){
        int tupleIndex1 = random(gen);
        int tupleIndex2 = random(gen);
        if (tupleIndex1 == tupleIndex2){
            i--;
            continue;
        }

        dynamic_bitset<> agreeSet(relationData->getNumColumns());
        for (auto & columnData : relationData->getColumnData()){
            int value1 = columnData->getProbingTableValue(tupleIndex1);
            if (value1 != PositionListIndex::singletonValueId && value1 == columnData->getProbingTableValue(tupleIndex2)){
                agreeSet[columnData->getColumn()->getIndex()] = true;
            }
        }

        agreeSetCounters[agreeSet]++;
    }

    T instance = std::make_shared<T>(relationData, relationData->getSchema()->emptyVertical, (int) sampleSize, relationData->getNumTuplePairs(), agreeSetCounters);
    return instance;
}

//template<typename T, typename enable_if<is_base_of<AgreeSetSample, T>::value>::type*>
template<typename T>
shared_ptr<T> AgreeSetSample::createFocusedFor(shared_ptr<ColumnLayoutRelationData> relation,
                                               shared_ptr<Vertical> restrictionVertical,
                                               shared_ptr<PositionListIndex> restrictionPli, unsigned int sampleSize) {
    static_assert(std::is_base_of<AgreeSetSample, T>::value);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> random_double;

    dynamic_bitset<> freeColumnIndices(relation->getNumColumns());
    freeColumnIndices.set();
    //auto restrictionVerticalNotIndices = ~restrictionVertical->getColumnIndices();    // ВОТ ТУТ ПОМЕНЯЛ. БЫЛО freeColumnIndices &=...; freeColumnIndices ~= freeColumnIndices
    freeColumnIndices &= ~restrictionVertical->getColumnIndices();
    vector<shared_ptr<ColumnData>> relevantColumnData;
    for (int columnIndex = freeColumnIndices.find_first(); columnIndex < freeColumnIndices.size(); columnIndex = freeColumnIndices.find_next(columnIndex)){
        relevantColumnData.push_back(relation->getColumnData(columnIndex));
    }
    dynamic_bitset<> agreeSetPrototype(restrictionVertical->getColumnIndices());
    map<dynamic_bitset<>, int> agreeSetCounters;

    unsigned long long restrictionNep = restrictionPli->getNepAsLong();
    sampleSize = std::min(static_cast<unsigned long long>(sampleSize), restrictionNep);
    if (sampleSize >= restrictionNep){
        for (auto & cluster : restrictionPli->getIndex()){
            for (int i = 0; i < cluster.size(); i++){
                int tupleIndex1 = cluster[i];
                for (int j = i + 1; j < cluster.size(); j++){
                    int tupleIndex2 = cluster[j];

                    dynamic_bitset<> agreeSet(agreeSetPrototype);
                    for (auto & columnData : relevantColumnData){
                        int value1 = columnData->getProbingTableValue(tupleIndex1);
                        if (value1 != PositionListIndex::singletonValueId && value1 == columnData->getProbingTableValue(tupleIndex2)){
                            agreeSet.set(columnData->getColumn()->getIndex());
                        }
                    }
                    if (agreeSetCounters.find(agreeSet) == agreeSetCounters.end()) {          // тут изменил
                        agreeSetCounters[agreeSet] = 1;
                    } else {
                        agreeSetCounters[agreeSet]++;
                    }
                }
            }
        }
    } else {
        vector<unsigned long long> clusterSizes(restrictionPli->getNumNonSingletonCluster() - 1);
        for (int i = 0; i < clusterSizes.size(); i++){
            unsigned long long clusterSize = restrictionPli->getIndex()[i].size();
            unsigned long long numTuplePairs = clusterSize * (clusterSize - 1) / 2;
            if (i > 0){
                clusterSizes[i] = numTuplePairs + clusterSizes[i - 1];
            } else {
                clusterSizes[i] = numTuplePairs;
            }
        }

        for (int i = 0; i < sampleSize; i++){
            auto clusterIndexIter = std::lower_bound(clusterSizes.begin(), clusterSizes.end(), (unsigned long long) (restrictionNep * random_double(gen)));
            unsigned int clusterIndex = std::distance(clusterSizes.begin(), clusterIndexIter);
            /*if (clusterIndex >= clusterSizes.size()) {
                clusterIndex = clusterSizes.size() - 1;
            }*/
            auto& cluster = restrictionPli->getIndex()[clusterIndex];
            std::uniform_int_distribution<> random_int(0, cluster.size() - 1);

            int tupleIndex1 = random_int(gen);
            int tupleIndex2 = random_int(gen);
            while (tupleIndex1 == tupleIndex2) {
                tupleIndex2 = random_int(gen);
            }
            tupleIndex1 = cluster[tupleIndex1];
            tupleIndex2 = cluster[tupleIndex2];


            dynamic_bitset<> agreeSet(agreeSetPrototype);
            for (auto & columnData : relevantColumnData){
                int value1 = columnData->getProbingTableValue(tupleIndex1);
                if (value1 != PositionListIndex::singletonValueId && value1 == columnData->getProbingTableValue(tupleIndex2)){
                    agreeSet.set(columnData->getColumn()->getIndex());
                }
            }

            if (agreeSetCounters.find(agreeSet) == agreeSetCounters.end()) {          // тут изменил
                agreeSetCounters[agreeSet] = 1;
            } else {
                agreeSetCounters[agreeSet]++;
            }
        }
    }
    /*std::cout << "-----------------\n";
    for (auto& [key, value] : agreeSetCounters) {
        //std::string tmp;
        //boost::to_string(key, tmp);
        //std::cout << '{';
        for (int columnIndex = key.find_first(); columnIndex < key.size(); columnIndex = key.find_next(columnIndex)){
            std::cout << columnIndex << ' ';
        }
        std::cout << '}';
        std::cout << '-' << value << '\n';
    }*/

    shared_ptr<T> sample = std::make_shared<T>(relation, restrictionVertical, sampleSize, restrictionNep, agreeSetCounters);
    return sample;
}
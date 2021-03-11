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

#include "ColumnLayoutRelationData.h"
#include "PositionListIndex.h"
#include "Vertical.h"


const int PositionListIndex::singletonValueId = 0;
unsigned long long PositionListIndex::micros = 0;
int PositionListIndex::intersectionCount = 0;

PositionListIndex::PositionListIndex(std::deque<std::vector<int>> index, std::vector<int> nullCluster,
                                     unsigned int size, double entropy, unsigned long long nep,
                                     unsigned int relationSize, unsigned int originalRelationSize,
                                     double invertedEntropy, double giniImpurity):
                                     index(std::move(index)),
                                     nullCluster(std::move(nullCluster)),
                                     size(size),
                                     entropy(entropy),
                                     invertedEntropy(invertedEntropy),
                                     giniImpurity(giniImpurity),
                                     nep(nep),
                                     relationSize(relationSize),
                                     originalRelationSize(originalRelationSize),
                                     probingTableCache(){}

std::unique_ptr<PositionListIndex> PositionListIndex::createFor(std::vector<int>& data, bool isNullEqNull) {
    std::unordered_map<int, std::vector<int>> index;
    for (unsigned long position = 0; position < data.size(); ++position){
        int valueId = data[position];
        index[valueId].push_back(position);
    }

    std::vector<int> nullCluster;
    if (index.count(RelationData::nullValueId) != 0) {
        nullCluster = index[RelationData::nullValueId];
    }
    if (!isNullEqNull){
        index.erase(RelationData::nullValueId); // move?
    }

    double keyGap = 0.0;
    double invEnt = 0;
    double giniGap = 0;
    unsigned long long nep = 0;
    unsigned int size = 0;
    std::deque<std::vector<int>> clusters;

    for (auto & iter : index){
        if (iter.second.size() == 1){
            giniGap += std::pow(1 / static_cast<double>(data.size()), 2);
            continue;
        }
        keyGap += iter.second.size() * log(iter.second.size());
        nep += calculateNep(iter.second.size());
        size += iter.second.size();
        invEnt += -(1 - iter.second.size() / static_cast<double>(data.size()))
                * std::log(1 - (iter.second.size() / static_cast<double>(data.size())));
        giniGap += std::pow(iter.second.size() / static_cast<double>(data.size()), 2);

        clusters.emplace_back(std::move(iter.second));
    }
    double entropy = log(data.size()) - keyGap / data.size();

    double giniImpurity = 1 - giniGap;
    if (giniImpurity == 0) {
        invEnt = 0;
    }

    sortClusters(clusters);
    return std::make_unique<PositionListIndex>(std::move(clusters), std::move(nullCluster),
                                               size, entropy, nep, data.size(), data.size(), invEnt, giniImpurity);
}

unsigned long long PositionListIndex::calculateNep(unsigned int numElements) {
    return static_cast<unsigned long long>(numElements) * (numElements - 1) / 2;
}

void PositionListIndex::sortClusters(std::deque<std::vector<int>> &clusters) {
    sort(clusters.begin(), clusters.end(), [](std::vector<int> & a, std::vector<int> & b){
        return a[0] < b[0]; } );
}

std::shared_ptr<const std::vector<int>> PositionListIndex::getProbingTable() const {
    if (probingTableCache != nullptr) return probingTableCache;

    std::vector<int> probingTable = std::vector<int>(originalRelationSize);
    int nextClusterId = singletonValueId + 1;
    for (auto & cluster : index){
        int valueId = nextClusterId++;
        assert(valueId != singletonValueId);
        for(int position : cluster){
            probingTable[position] = valueId;
        }
    }

    return std::make_shared<std::vector<int>>(probingTable);
}

void PositionListIndex::forceCacheProbingTable() {
    getProbingTable(true);
}

// интересное место: true --> надо передать поле без копирования, false --> надо сконструировать и выдать наружу
// кажется, самым лёгким способом будет навернуть shared_ptr
std::shared_ptr<const std::vector<int>> PositionListIndex::getProbingTable(bool isCaching) {
    auto probingTable = getProbingTable();
    if (isCaching) {
        probingTableCache = probingTable;
        return probingTableCache;
    }
    return probingTable;
}

std::deque<std::vector<int>> const & PositionListIndex::getIndex() const {
    return index;
}

std::unique_ptr<PositionListIndex> PositionListIndex::intersect(PositionListIndex const* that) const {
    assert(this->relationSize == that->relationSize);
    return this->size > that->size ?
            that->probe(this->getProbingTable()) :
            this->probe(that->getProbingTable());
}

//TODO: nullCluster некорректен
std::unique_ptr<PositionListIndex> PositionListIndex::probe(std::shared_ptr<const std::vector<int>> probingTable) const {
    assert(this->relationSize == probingTable->size());
    std::deque<std::vector<int>> newIndex;
    unsigned int newSize = 0;
    double newKeyGap = 0.0;
    unsigned long long newNep = 0;
    std::vector<int> nullCluster;

    std::unordered_map<int, std::vector<int>> partialIndex;

    for (auto & positions : index){
        for (int position : positions){
            int probingTableValueId = (*probingTable)[position];
            if (probingTableValueId == singletonValueId)
                continue;
            intersectionCount++;
            partialIndex[probingTableValueId].push_back(position);
        }

        for (auto & iter : partialIndex){
            auto & cluster = iter.second;
            if (cluster.size() <= 1) continue;

            newSize += cluster.size();
            newKeyGap += cluster.size() * log(cluster.size());
            newNep += calculateNep(cluster.size());

            newIndex.push_back(std::move(cluster));
        }
        partialIndex.clear();
    }

    double newEntropy = log(relationSize) - newKeyGap / relationSize;
    sortClusters(newIndex);

    return std::make_unique<PositionListIndex>(std::move(newIndex), std::move(nullCluster),
                                               newSize, newEntropy, newNep, relationSize, relationSize);
}


//TODO: nullCluster не поддерживается
std::unique_ptr<PositionListIndex> PositionListIndex::probeAll(Vertical const& probingColumns,
                                                               ColumnLayoutRelationData & relationData) {
    assert(this->relationSize == relationData.getNumRows());
    std::deque<std::vector<int>> newIndex;
    unsigned int newSize = 0;
    double newKeyGap = 0.0;
    unsigned long long newNep = 0;

    std::map<std::vector<int>, std::vector<int>> partialIndex;
    std::vector<int> nullCluster;
    std::vector<int> probe;

    for (auto & cluster : this->index){
        for (int position : cluster){
            if (!takeProbe(position, relationData, probingColumns, probe)){
                probe.clear();
                continue;
            }

            partialIndex[probe].push_back(position);
            probe.clear();
        }

        for (auto & iter : partialIndex) {
            auto & newCluster = iter.second;
            if (newCluster.size() == 1) continue;

            newSize += newCluster.size();
            newKeyGap += newCluster.size() * log(newCluster.size());
            newNep += calculateNep(newCluster.size());

            newIndex.emplace_back(std::move(newCluster));
        }
        partialIndex.clear();
    }

    double newEntropy = log(this->relationSize) - newKeyGap / this->relationSize;

    sortClusters(newIndex);

    return std::make_unique<PositionListIndex>(
            std::move(newIndex), std::move(nullCluster), newSize, newEntropy, newNep, this->relationSize, this->relationSize
            );
}

bool PositionListIndex::takeProbe(int position, ColumnLayoutRelationData & relationData,
                                  Vertical const& probingColumns, std::vector<int> & probe) {
    boost::dynamic_bitset<> probingIndices = probingColumns.getColumnIndices();
    for (unsigned long index = probingIndices.find_first();
         index < probingIndices.size();
         index = probingIndices.find_next(index)){
        int value = relationData.getColumnData(index).getProbingTableValue(position);
        if (value == PositionListIndex::singletonValueId) return false;
        probe.push_back(value);
    }
    return true;
}

std::string PositionListIndex::toString() const {
    std::string res = "[";
    for (auto& cluster : index){
        res.push_back('[');
        for (int v : cluster){
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



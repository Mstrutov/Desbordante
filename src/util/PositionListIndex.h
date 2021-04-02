//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <memory>
#include <deque>
#include <vector>

#include "Column.h"

using std::deque, std::vector, std::shared_ptr;

class ColumnLayoutRelationData;

class PositionListIndex {
private:
    deque<vector<int>> index;
    vector<int> nullCluster;
    unsigned int size;
    double entropy;
    double invertedEntropy;
    double giniImpurity;
    unsigned long long nep;
    unsigned int relationSize;
    unsigned int originalRelationSize;
    vector<int> probingTableCache;
    unsigned int freq_ = 0;
    int clusterNum = -1;

    static unsigned long long calculateNep(unsigned int numElements);
    static void sortClusters(deque<vector<int>> & clusters);
    static bool takeProbe(int position, ColumnLayoutRelationData & relationData, Vertical & probingColumns, vector<int> & probe);

public:
    PositionListIndex(deque<vector<int>>&& index, vector<int>&& nullCluster, unsigned int size, double entropy,
                      unsigned long long nep, unsigned int relationSize, unsigned int originalRelationSize,
                      double invertedEntropy = 0, double giniImpurity = 0);
    //TODO: for optimization needs
    static int intersectionCount;
    static unsigned long long micros;
    static const int singletonValueId;
    static shared_ptr<PositionListIndex> createFor(vector<int>& data, bool isNullEqNull);
    vector<int> getProbingTable();
    vector<int> getProbingTable(bool isCaching);
    deque<vector<int>> const & getIndex();
    double getNep() { return (double) nep; }
    unsigned long long getNepAsLong() const { return nep; }           //ADDED: getNep()
    unsigned int getNumNonSingletonCluster() const { return index.size(); }
    unsigned int getNumCluster();
    unsigned int getFreq() const { return freq_; }
    void incFreq() { freq_++; }
    int getSize() { return size; }
    double getEntropy() const { return entropy; }
    double getInvertedEntropy() const { return invertedEntropy; }
    double getGiniImpurity() const { return giniImpurity; }
    double getMaximumNip() const { return calculateNep(relationSize); }
    double getNip() const { return getMaximumNip() - getNepAsLong(); }

    shared_ptr<PositionListIndex> intersect(shared_ptr<PositionListIndex> that);
    shared_ptr<PositionListIndex> probe(const vector<int>& probingTable);
    shared_ptr<PositionListIndex> probeAll(Vertical probingColumns, ColumnLayoutRelationData & relationData);

    std::string toString() const;
};

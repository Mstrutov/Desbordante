//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <memory>
#include <deque>
#include <vector>

#include "Column.h"

class ColumnLayoutRelationData;

class PositionListIndex {
private:
    std::deque<std::vector<int>> index;
    std::vector<int> nullCluster;
    unsigned int size;
    double entropy;
    double invertedEntropy;
    double giniImpurity;
    unsigned long long nep;
    unsigned int relationSize;
    unsigned int originalRelationSize;
    std::shared_ptr<const std::vector<int>> probingTableCache;
    unsigned int freq_ = 0;


    static unsigned long long calculateNep(unsigned int numElements);
    static void sortClusters(std::deque<std::vector<int>> & clusters);
    static bool takeProbe(int position, ColumnLayoutRelationData & relationData,
                          Vertical const& probingColumns, std::vector<int> & probe);

public:
    static int intersectionCount;
    static unsigned long long micros;
    static const int singletonValueId;

    PositionListIndex(std::deque<std::vector<int>> index, std::vector<int> nullCluster,
                      unsigned int size, double entropy,
                      unsigned long long nep, unsigned int relationSize, unsigned int originalRelationSize,
                      double invertedEntropy = 0, double giniImpurity = 0);
    static std::unique_ptr<PositionListIndex> createFor(std::vector<int>& data, bool isNullEqNull);

    std::shared_ptr<const std::vector<int>> getProbingTable() const;
    void forceCacheProbingTable();
    std::shared_ptr<const std::vector<int>> getProbingTable(bool isCaching);

    std::deque<std::vector<int>> const & getIndex() const;
    double getNep()                             const { return (double) nep; }
    unsigned long long getNepAsLong()           const { return nep; }
    unsigned int getNumNonSingletonCluster()    const { return index.size(); }
    unsigned int getFreq()                      const { return freq_; }
    unsigned int getSize()                      const { return size; }
    double getEntropy()                         const { return entropy; }
    double getInvertedEntropy()                 const { return invertedEntropy; }
    double getGiniImpurity()                    const { return giniImpurity; }
    double getMaximumNip()                      const { return calculateNep(relationSize); }
    double getNip()                             const { return getMaximumNip() - getNepAsLong(); }

    void incFreq() { freq_++; }

    std::unique_ptr<PositionListIndex> intersect(PositionListIndex const* that) const;
    std::unique_ptr<PositionListIndex> probe(std::shared_ptr<const std::vector<int>> probingTable) const;
    std::unique_ptr<PositionListIndex> probeAll(Vertical const& probingColumns,
                                                ColumnLayoutRelationData & relationData);
    std::string toString() const;
};

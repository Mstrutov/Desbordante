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


    static unsigned long long calculateNep(unsigned int numElements) {
        return static_cast<unsigned long long>(numElements) * (numElements - 1) / 2;
    }
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

    // если PT закеширована, выдаёт её, иначе предварительно вычисляет её -- тяжёлая операция
    std::shared_ptr<const std::vector<int>> calculateAndGetProbingTable() const;
    // выдаёт закешированную PT, либо nullptr, если она не закеширована
    std::vector<int> const* getCachedProbingTable() const { return probingTableCache.get(); };
    // кеширует PT
    void forceCacheProbingTable() { probingTableCache = calculateAndGetProbingTable(); };
    // Такая структура с кешированием ProbingTable нужна, потому что к PT одиночных колонок происходят
    // частые обращения, чтобы узнать какую-то одну конкретную позицию, тогда как PT наборов колонок
    // обычно используются, чтобы один раз пересечь две партиции, и больше к ним не возвращаться

    // std::shared_ptr<const std::vector<int>> getProbingTable(bool isCaching);

    std::deque<std::vector<int>> const & getIndex() const { return index; };
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

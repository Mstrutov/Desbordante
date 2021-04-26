#pragma once

#include <set>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "FDAlgorithm.h"
#include "CSVParser.h"
#include "ColumnCombination.h"
#include "ColumnLayoutRelationData.h"
#include "PositionListIndex.h"
#include "Vertical.h"

class Fd_mine : public FDAlgorithm {
  private:

    shared_ptr<ColumnLayoutRelationData> relation;
    shared_ptr<RelationalSchema> schema;

    std::set<dynamic_bitset<>> candidateSet;
    boost::unordered_map<dynamic_bitset<>, std::unordered_set<dynamic_bitset<>>> eqSet;
    boost::unordered_map<dynamic_bitset<>, dynamic_bitset<>> fdSet;
    std::set<dynamic_bitset<>> keySet;
    boost::unordered_map<dynamic_bitset<>, dynamic_bitset<>> closure;
    boost::unordered_map<dynamic_bitset<>, shared_ptr<PositionListIndex>> plis;
    dynamic_bitset<> r;

    void computeNonTrivialClosure(dynamic_bitset<> xi);
    void obtainFDandKey(dynamic_bitset<> xi);
    void obtainEQSet();
    void pruneCandidates();
    void generateCandidates();
    void display();
    
  public:
    Fd_mine(fs::path const& path, char separator = ',', bool hasHeader = true) : FDAlgorithm(path, separator, hasHeader){};
    ~Fd_mine() override {};
    unsigned long long execute() override;
};

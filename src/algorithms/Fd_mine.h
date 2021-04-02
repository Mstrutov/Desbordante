#pragma once

#include <set>

#include "FDAlgorithm.h"
#include "CSVParser.h"
#include "ColumnCombination.h"
#include "ColumnLayoutRelationData.h"
#include "PositionListIndex.h"

class Fd_mine : public FDAlgorithm {
  private:

    shared_ptr<ColumnLayoutRelationData> relation;
    shared_ptr<RelationalSchema> schema;

    std::set<dynamic_bitset<>> candidateSet;
    std::map<dynamic_bitset<>, std::set<dynamic_bitset<>>> eqSet;
    std::map<dynamic_bitset<>, dynamic_bitset<>> fdSet;
    std::set<dynamic_bitset<>> keySet;
    std::map<dynamic_bitset<>, dynamic_bitset<>> closure;
    std::map<dynamic_bitset<>, shared_ptr<PositionListIndex>> plis;
    dynamic_bitset<> r;

    void computeNonTrivialClosure(dynamic_bitset<> xi);
    void obtainFDandKey(dynamic_bitset<> xi);
    void obtainEQSet(std::set<dynamic_bitset<>> &candidateSet);
    void pruneCandidates(std::set<dynamic_bitset<>> &candidateSet);
    void generateCandidates(std::set<dynamic_bitset<>> &candidateSet);
    void display();
    
  public:
    Fd_mine(fs::path& path) : FDAlgorithm(path){};
    ~Fd_mine() override {};
    unsigned long long execute();
};

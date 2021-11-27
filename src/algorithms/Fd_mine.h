#pragma once

#include <boost/unordered_map.hpp>
#include <filesystem>
#include <set>

#include "CSVParser.h"
#include "ColumnCombination.h"
#include "ColumnLayoutRelationData.h"
#include "PliBasedFDAlgorithm.h"
#include "PositionListIndex.h"
#include "Vertical.h"

class Fd_mine : public PliBasedFDAlgorithm {
   private:
    const RelationalSchema* schema;

    std::set<dynamic_bitset<>> candidateSet;
    boost::unordered_map<dynamic_bitset<>, std::unordered_set<dynamic_bitset<>>> eqSet;
    boost::unordered_map<dynamic_bitset<>, dynamic_bitset<>> fdSet;
    boost::unordered_map<dynamic_bitset<>, dynamic_bitset<>> final_fdSet;
    std::set<dynamic_bitset<>> keySet;
    boost::unordered_map<dynamic_bitset<>, dynamic_bitset<>> closure;
    boost::unordered_map<dynamic_bitset<>, std::shared_ptr<util::PositionListIndex const>> plis;
    dynamic_bitset<> relationIndices;

    void computeNonTrivialClosure(dynamic_bitset<> const& xi);
    void obtainFDandKey(dynamic_bitset<> const& xi);
    void obtainEQSet();
    void pruneCandidates();
    void generateNextLevelCandidates();
    void reconstruct();
    void display();

    unsigned long long executeInternal() override;
public:
    Fd_mine(std::filesystem::path const& path, char separator = ',', bool hasHeader = true)
            : PliBasedFDAlgorithm(path, separator, hasHeader) {}
    ~Fd_mine() override {}
};

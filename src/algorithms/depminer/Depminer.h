#pragma once

#include "CSVParser.h"
#include "FDAlgorithm.h"
#include "PliBasedFDAlgorithm.h"
#include "depminer/CMAXSet.h"

class Depminer : public PliBasedFDAlgorithm {
private:
    static CMAXSet GenFirstLevel(std::vector<CMAXSet> const& cmax_sets, Column const& attribute,
                                 std::unordered_set<Vertical>& level);
    static std::unordered_set<Vertical> GenNextLevel(
        std::unordered_set<Vertical> const& prev_level);
    static bool CheckJoin(Vertical const& _p, Vertical const& _q);

    void LhsForColumn(std::unique_ptr<Column> const& column, std::vector<CMAXSet> const& cmax_sets);
    std::vector<CMAXSet> GenerateCmaxSets(std::unordered_set<Vertical> const& agree_sets);

    double progress_step_ = 0;

public:
    explicit Depminer(std::filesystem::path const& path, char separator = ',',
                      bool has_header = true)
        : PliBasedFDAlgorithm(path, separator, has_header, true,
                              {"AgreeSets generation", "Finding CMAXSets", "Finding LHS"}) {}
    explicit Depminer(std::shared_ptr<ColumnLayoutRelationData> relation)
        : PliBasedFDAlgorithm(std::move(relation),
                              {"AgreeSets generation", "Finding CMAXSets", "Finding LHS"}) {}

    unsigned long long ExecuteInternal() override;

    const RelationalSchema* schema_ = nullptr;
};

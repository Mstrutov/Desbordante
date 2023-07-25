#pragma once

#include "algorithms/depminer/cmax_set.h"
#include "algorithms/pli_based_fd_algorithm.h"

namespace algos {

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
    RelationalSchema const* schema_ = nullptr;

    void ResetStateFd() final {}
    unsigned long long ExecuteInternal() final;

public:
    Depminer(bool request_prepared_data = false);
};

}  // namespace algos

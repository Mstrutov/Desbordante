#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "algorithms/md/hymd/indexes/records_info.h"
#include "algorithms/md/hymd/lattice/full_lattice.h"
#include "algorithms/md/hymd/preprocessing/similarity_measure/similarity_measure.h"
#include "algorithms/md/hymd/similarity_data.h"
#include "algorithms/md/hymd/similarity_measure_creator.h"
#include "algorithms/md/md_algorithm.h"
#include "config/tabular_data/input_table_type.h"
#include "model/table/relational_schema.h"

namespace algos::hymd {

class HyMD final : public MdAlgorithm {
private:
    using ColMatchesVector = std::vector<
            std::tuple<std::string, std::string, std::shared_ptr<SimilarityMeasureCreator>>>;

    config::InputTable left_table_;
    config::InputTable right_table_;

    std::shared_ptr<RelationalSchema> left_schema_;
    std::shared_ptr<RelationalSchema> right_schema_;

    std::unique_ptr<indexes::RecordsInfo> records_info_;

    std::size_t min_support_ = 0;
    bool prune_nondisjoint_ = true;
    // TODO: thread number limit
    // TODO: different level definitions (cardinality currently used)
    // TODO: comparing only some values during similarity calculation
    // TODO: cardinality limit
    // TODO: automatically calculating minimal support
    // TODO: limit LHS bounds searched (currently only size limit is implemented)
    // TODO: memory conservation mode (load only some columns)

    ColMatchesVector column_matches_option_;

    void RegisterOptions();

    void LoadDataInternal() final;

    void MakeExecuteOptsAvailable() final;
    void ResetStateMd() final;
    unsigned long long ExecuteInternal() final;

    void RegisterResults(SimilarityData const& similarity_data,
                         std::vector<lattice::MdLatticeNodeInfo> lattice_mds);

public:
    HyMD();
};

}  // namespace algos::hymd

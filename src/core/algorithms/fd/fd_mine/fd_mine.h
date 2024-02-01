#pragma once

#include <filesystem>
#include <set>

#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>

#include "algorithms/fd/pli_based_fd_algorithm.h"
#include "model/table/column_layout_relation_data.h"
#include "model/table/position_list_index.h"
#include "model/table/vertical.h"

namespace algos {

class Fd_mine : public PliBasedFDAlgorithm {
private:
    RelationalSchema const* schema_;

    std::set<boost::dynamic_bitset<>> candidate_set_;
    boost::unordered_map<boost::dynamic_bitset<>, std::unordered_set<boost::dynamic_bitset<>>>
            eq_set_;
    boost::unordered_map<boost::dynamic_bitset<>, boost::dynamic_bitset<>> fd_set_;
    boost::unordered_map<boost::dynamic_bitset<>, boost::dynamic_bitset<>> final_fd_set_;
    std::set<boost::dynamic_bitset<>> key_set_;
    boost::unordered_map<boost::dynamic_bitset<>, boost::dynamic_bitset<>> closure_;
    boost::unordered_map<boost::dynamic_bitset<>, std::shared_ptr<model::PositionListIndex const>>
            plis_;
    boost::dynamic_bitset<> relation_indices_;

    void ComputeNonTrivialClosure(boost::dynamic_bitset<> const& xi);
    void ObtainFDandKey(boost::dynamic_bitset<> const& xi);
    void ObtainEqSet();
    void PruneCandidates();
    void GenerateNextLevelCandidates();
    void Reconstruct();
    void Display();

    void ResetStateFd() final;
    unsigned long long ExecuteInternal() override;

public:
    Fd_mine(std::optional<ColumnLayoutRelationDataManager> relation_manager = std::nullopt);
};

}  // namespace algos

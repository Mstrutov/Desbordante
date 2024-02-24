#pragma once

#include <string>
#include <vector>

#include "algorithms/fd/fd_algorithm.h"
#include "algorithms/fd/fdep/fd_tree_element.h"
#include "config/equal_nulls/type.h"
#include "config/tabular_data/input_table_type.h"
#include "model/table/relation_data.h"
#include "model/table/relational_schema.h"

namespace algos {

class FDep : public FDAlgorithm {
public:
    FDep();

    ~FDep() override = default;

private:
    config::InputTable input_table_;

    std::unique_ptr<RelationalSchema> schema_{};

    std::vector<std::string> column_names_;
    size_t number_attributes_{};

    std::unique_ptr<FDTreeElement> neg_cover_tree_{};
    std::unique_ptr<FDTreeElement> pos_cover_tree_{};

    std::vector<std::vector<size_t>> tuples_;

    void RegisterOptions();

    void LoadDataInternal() final;

    void ResetStateFd() final;
    unsigned long long ExecuteInternal() final;

    // Building negative cover via violated dependencies
    void BuildNegativeCover();

    // Iterating over all pairs t1 and t2 of the relation
    // Adding violated FDs to negative cover tree.
    void AddViolatedFDs(std::vector<size_t> const& t1, std::vector<size_t> const& t2);

    // Converting negative cover tree into positive cover tree
    void CalculatePositiveCover(FDTreeElement const& neg_cover_subtree,
                                std::bitset<FDTreeElement::kMaxAttrNum>& active_path);

    // Specializing general dependencies for not to be followed from violated dependencies of
    // negative cover tree.
    void SpecializePositiveCover(std::bitset<FDTreeElement::kMaxAttrNum> const& lhs,
                                 size_t const& a);
};

}  // namespace algos

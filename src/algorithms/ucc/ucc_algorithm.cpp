#include "algorithms/ucc/ucc_algorithm.h"

#include "util/config/equal_nulls/option.h"
#include "util/config/tabular_data/input_table/option.h"

namespace algos {

UCCAlgorithm::UCCAlgorithm(std::vector<std::string_view> phase_names)
    : Algorithm(std::move(phase_names)) {
    RegisterOptions();
}

void UCCAlgorithm::RegisterOptions() {
    RegisterInitialLoadOption(util::config::TableOpt(&input_table_));
    RegisterInitialLoadOption(util::config::EqualNullsOpt(&is_null_equal_null_));
}

}  // namespace algos

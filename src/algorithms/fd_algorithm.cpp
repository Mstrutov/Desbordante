#include "algorithms/fd_algorithm.h"

#include <thread>

#include "util/config/equal_nulls/option.h"
#include "util/config/tabular_data/input_table/option.h"

namespace algos {

FDAlgorithm::FDAlgorithm(std::vector<std::string_view> phase_names,
                         util::config::ConfigurationStage initial_stage)
    : Algorithm(std::move(phase_names), initial_stage) {
    RegisterOptions();
}

void FDAlgorithm::RegisterOptions() {
    RegisterInitialLoadOption(util::config::TableOpt(&input_table_));
    RegisterInitialLoadOption(util::config::EqualNullsOpt(&is_null_equal_null_));
}

void FDAlgorithm::ResetState() {
    fd_collection_.Clear();
    ResetStateFd();
}

std::string FDAlgorithm::GetJsonFDs() const {
    return FDsToJson(FdList());
}

unsigned int FDAlgorithm::Fletcher16() {
    std::string to_hash = GetJsonFDs();
    unsigned int sum1 = 0, sum2 = 0, modulus = 255;
    for (auto ch : to_hash) {
        sum1 = (sum1 + ch) % modulus;
        sum2 = (sum2 + sum1) % modulus;
    }
    return (sum2 << 8) | sum1;
}

/* Attribute A contains only unique values (i.e. A is the key) iff [A]->[B]
 * holds for every attribute B. So to determine if A is a key, we count
 * number of fds with lhs==[A] and if it equals the total number of attributes
 * minus one (the attribute A itself) then A is the key.
 */
std::vector<Column const*> FDAlgorithm::GetKeys() const {
    std::vector<Column const*> keys;
    std::map<Column const*, size_t> fds_count_per_col;
    unsigned int cols_of_equal_values = 0;

    for (FD const& fd : FdList()) {
        Vertical const& lhs = fd.GetLhs();

        if (lhs.GetArity() == 0) {
            /* We separately count columns consisting of only equal values,
             * because they cannot be on the right side of the minimal fd.
             * And obviously for every attribute A true: [A]->[B] holds
             * if []->[B] holds.
             */
            cols_of_equal_values++;
        } else if (lhs.GetArity() == 1) {
            fds_count_per_col[lhs.GetColumns().front()]++;
        }
    }

    if (fds_count_per_col.empty()) return keys;
    assert(fds_count_per_col.begin()->first->GetSchema() != nullptr);
    size_t const number_of_columns = fds_count_per_col.begin()->first->GetSchema()->GetNumColumns();
    [[maybe_unused]] RelationalSchema const* first_schema =
            fds_count_per_col.begin()->first->GetSchema();
    for (auto const& [col, num]: fds_count_per_col) {
        assert(col->GetSchema() == first_schema);
        if (num + 1 + cols_of_equal_values == number_of_columns) {
            keys.push_back(col);
        }
    }

    return keys;
}

}  // namespace algos

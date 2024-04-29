#include "fd_algorithm.h"

#include <map>
#include <thread>
#include <vector>

#include "config/max_lhs/option.h"
#include "model/table/relational_schema.h"

namespace algos {

FDAlgorithm::FDAlgorithm(std::vector<std::string_view> phase_names)
    : Algorithm(std::move(phase_names)) {
    RegisterOptions();
}

void FDAlgorithm::RegisterOptions() {
    RegisterOption(config::kMaxLhsOpt(&max_lhs_));
}

void FDAlgorithm::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable({config::kMaxLhsOpt.GetName()});
    MakeExecuteOptsAvailableFDInternal();
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
    for (auto const& [col, num] : fds_count_per_col) {
        assert(col->GetSchema() == first_schema);
        if (num + 1 + cols_of_equal_values == number_of_columns) {
            keys.push_back(col);
        }
    }

    return keys;
}

}  // namespace algos

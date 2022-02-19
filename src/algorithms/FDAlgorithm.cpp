#include "FDAlgorithm.h"

#include <thread>

unsigned long long FDAlgorithm::Execute() {
    Initialize();

    return ExecuteInternal();
}

void FDAlgorithm::InitConfigParallelism() {
    if (config_.parallelism == 0) {
        config_.parallelism = std::thread::hardware_concurrency();
        if (config_.parallelism == 0) {
            throw std::runtime_error("Unable to detect number of concurrent "
                                     "threads supported by your system. "
                                     "Please, specify it manually.");
        }
    }
}

std::string FDAlgorithm::GetJsonFDs() const {
        return FDsToJson(fd_collection_);
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
    size_t const number_of_cols = input_generator_.GetNumberOfColumns();

    for (FD const& fd : fd_collection_) {
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

    for (auto const&[col, num] : fds_count_per_col) {
        if (num + 1 + cols_of_equal_values == number_of_cols) {
            keys.push_back(col);
        }
    }

    return keys;
}


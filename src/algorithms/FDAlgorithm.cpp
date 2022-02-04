#include "FDAlgorithm.h"

unsigned long long FDAlgorithm::Execute() {
    Initialize();

    return ExecuteInternal();
}

std::string FDAlgorithm::GetJsonFDs() {
    std::string result = "{\"fds\": [";
    std::list<std::string> discovered_fd_strings;
    for (auto& fd : fd_collection_) {
        discovered_fd_strings.push_back(fd.ToJSONString());
    }
    discovered_fd_strings.sort();
    for (auto const& fd : discovered_fd_strings) {
        result += fd + ",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }
    /*result += ']';

    result += ", \"uccs\": [";
    std::list<std::string> discoveredUCCStrings;
    for (auto& ucc : discovered_uccs_) {
        discoveredUCCStrings.push_back(ucc.ToIndicesString());
    }
    discoveredUCCStrings.sort();
    for (auto const& ucc : discoveredUCCStrings) {
        result += '\"' + ucc + "\",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }*/
    result += "]}";
    return result;
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

void FDAlgorithm::AddProgress(double const val) noexcept {
        assert(val >= 0);
        std::scoped_lock lock(progress_mutex_);
        cur_phase_progress_ += val;
        assert(cur_phase_progress_ < 101);
}

void FDAlgorithm::SetProgress(double const val) noexcept {
        assert(0 <= val && val < 101);
        std::scoped_lock lock(progress_mutex_);
        cur_phase_progress_ = val;
}

std::pair<uint8_t, double> FDAlgorithm::GetProgress() const noexcept {
        std::scoped_lock lock(progress_mutex_);
        return std::make_pair(cur_phase_id_, cur_phase_progress_);
}

void FDAlgorithm::ToNextProgressPhase() noexcept {
    std::scoped_lock lock(progress_mutex_);
    ++cur_phase_id_;
    assert(cur_phase_id_ < phase_names_.size());
    assert(cur_phase_progress_ >= 100);
    cur_phase_progress_ = 0;
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


#include "FDAlgorithm.h"

std::string FDAlgorithm::getJsonFDs() {
    std::string result = "{\"fds\": [";
    std::list<std::string> discoveredFDStrings;
    for (auto& fd : fdCollection_) {
        discoveredFDStrings.push_back(fd.toJSONString());
    }
    discoveredFDStrings.sort();
    for (auto const& fd : discoveredFDStrings) {
        result += fd + ",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }
    /*result += ']';

    result += ", \"uccs\": [";
    std::list<std::string> discoveredUCCStrings;
    for (auto& ucc : discoveredUCCs_) {
        discoveredUCCStrings.push_back(ucc.toIndicesString());
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

unsigned int FDAlgorithm::fletcher16() {
    std::string toHash = getJsonFDs();
    unsigned int sum1 = 0, sum2 = 0, modulus = 255;
    for (auto ch : toHash) {
        sum1 = (sum1 + ch) % modulus;
        sum2 = (sum2 + sum1) % modulus;
    }
    return (sum2 << 8) | sum1;
}

void FDAlgorithm::addProgress(double const val) noexcept {
        assert(val >= 0);
#ifdef __cpp_lib_atomic_float
        progress_.fetch_add(val);
        assert(progress_.load() < 101);
#else
        std::scoped_lock lock(progress_mutex_);
        progress_ += val;
        assert(progress_ < 101);
#endif
}

void FDAlgorithm::setProgress(double const val) noexcept {
        assert(0 <= val && val < 101);
#ifdef __cpp_lib_atomic_float
        progress_.store(val);
#else
        std::scoped_lock lock(progress_mutex_);
        progress_ = val;
#endif
}

double FDAlgorithm::getProgress() const noexcept {
#ifdef __cpp_lib_atomic_float
        return progress_.load();
#else
        std::scoped_lock lock(progress_mutex_);
        return progress_;
#endif
}

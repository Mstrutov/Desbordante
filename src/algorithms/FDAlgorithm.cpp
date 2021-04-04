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


#include "DependencyConsumer.h"

PartialFD DependencyConsumer::registerFd(Vertical const& lhs, Column const& rhs,
                                         double error, double score) const {
    PartialFD partialFd(lhs, rhs, error, score);
    fdConsumer_(partialFd);
    return partialFd;
}

PartialKey DependencyConsumer::registerUcc(Vertical const& keyVertical,
                                           double error, double score) const {
    PartialKey partialKey(keyVertical, error, score);
    uccConsumer_(partialKey);
    return partialKey;
}

std::string DependencyConsumer::fdsToString() const {
    std::string result;
    for (auto const& fd : discoveredFDs_) {
        result += fd.toString() + "\n\r";
    }
    return result;
}

std::string DependencyConsumer::uccsToString() const {
    std::string result;
    for (auto const& ucc : discoveredUCCs_) {
        result += ucc.toString() + "\n\r";
    }
    return result;
}

std::string DependencyConsumer::getJsonFDs() {
    std::string result = "{\"fds\": [";
    std::list<std::string> discoveredFDStrings;
    for (auto& fd : discoveredFDs_) {
        discoveredFDStrings.push_back(fd.toIndicesString());
    }
    discoveredFDStrings.sort();
    for (auto const& fd : discoveredFDStrings) {
        result += '\"' + fd + "\",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }
    result += ']';

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
    }
    result += "]}";
    return result;
}


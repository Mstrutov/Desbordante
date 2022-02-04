
#include "DependencyConsumer.h"

PartialFD DependencyConsumer::RegisterFd(Vertical const& lhs, Column const& rhs,
                                         double error, double score) const {
    PartialFD partial_fd(lhs, rhs, error, score);
    fd_consumer_(partial_fd);
    return partial_fd;
}

PartialKey DependencyConsumer::RegisterUcc(Vertical const& key_vertical,
                                           double error, double score) const {
    PartialKey partial_key(key_vertical, error, score);
    ucc_consumer_(partial_key);
    return partial_key;
}

std::string DependencyConsumer::FDsToString() const {
    std::string result;
    for (auto const& fd : discovered_fds_) {
        result += fd.ToString() + "\n\r";
    }
    return result;
}

std::string DependencyConsumer::UCCsToString() const {
    std::string result;
    for (auto const& ucc : discovered_uccs_) {
        result += ucc.ToString() + "\n\r";
    }
    return result;
}

std::string DependencyConsumer::GetJsonFDs() {
    std::string result = "{\"fds\": [";
    std::list<std::string> discovered_fd_strings;
    for (auto& fd : discovered_fds_) {
        discovered_fd_strings.push_back(fd.ToIndicesString());
    }
    discovered_fd_strings.sort();
    for (auto const& fd : discovered_fd_strings) {
        result += '\"' + fd + "\",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }
    result += ']';

    result += ", \"uccs\": [";
    std::list<std::string> discovered_ucc_strings;
    for (auto& ucc : discovered_uccs_) {
        discovered_ucc_strings.push_back(ucc.ToIndicesString());
    }
    discovered_ucc_strings.sort();
    for (auto const& ucc : discovered_ucc_strings) {
        result += '\"' + ucc + "\",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }
    result += "]}";
    return result;
}


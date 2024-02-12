#include "raw_cfd.h"

#include <sstream>

namespace algos::cfd {

static std::string RawItemToJSON(RawCFD::RawItem const& raw_item) {
    std::stringstream ss;
    ss << "{\"attribute\":" << raw_item.attribute;
    ss << ",\"value\":";
    if (raw_item.value.has_value()) {
        ss << '"' << raw_item.value.value() << '"';
    } else {
        ss << "null";
    }
    ss << "}";
    return ss.str();
}

std::string RawCFD::ToJSON() const {
    std::stringstream ss;
    ss << "{\"lhs\":";
    ss << "[";
    for (auto it = lhs_.begin(); it != lhs_.end(); ++it) {
        if (it != lhs_.begin()) {
            ss << ",";
        }
        ss << RawItemToJSON(*it);
    }
    ss << "]";
    ss << ",\"rhs\":" + RawItemToJSON(rhs_) + "}";
    return ss.str();
}

static std::string RawItemToString(RawCFD::RawItem const& item) {
    std::stringstream ss;
    ss << "(" << item.attribute;
    ss << ", ";
    if (item.value.has_value()) {
        ss << item.value.value();
    } else {
        ss << "_";
    }
    ss << ")";
    return ss.str();
}

std::string RawCFD::ToString() const {
    std::stringstream ss;
    ss << "{";
    for (auto it = lhs_.begin(); it != lhs_.end(); ++it) {
        if (it != lhs_.begin()) {
            ss << ",";
        }
        ss << RawItemToString(*it);
    }
    ss << "} -> " << RawItemToString(rhs_);
    return ss.str();
}

}  // namespace algos::cfd

#pragma once

#include <sstream>
#include <string>
#include <vector>

namespace algos::nd_verifier::util {

template <typename T>
concept Printable = requires(std::stringstream& sstream, T& t) {
    sstream << t;
};

template <Printable T>
std::string VectorToString(std::vector<T> const& vect) {
    std::stringstream sstream;
    sstream << '[';
    for (auto pt{vect.begin()}; pt != vect.end(); ++pt) {
        if (pt != vect.begin()) {
            sstream << ", ";
        }
        sstream << *pt;
    }
    sstream << ']';
    return sstream.str();
}

}  // namespace algos::nd_verifier::util

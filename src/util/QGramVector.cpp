#include "QGramVector.h"

#include <cassert>
#include <cmath>
#include <numeric>

namespace util {

QGramVector::QGramVector(std::string_view const& string, unsigned q) {
    assert(string.size() >= q);
    for (size_t i = 0; i < string.size() - q + 1; ++i) {
        q_grams_[std::string(string.substr(i, q))]++;
    }
    CalculateLength();
}

double QGramVector::InnerProduct(QGramVector const& other) const {
    auto const& q_grams_pair = std::minmax(
        q_grams_, other.q_grams_, [](auto const& a, auto const& b) { return a.size() < b.size(); });
    auto const& l = q_grams_pair.first;
    auto const& r = q_grams_pair.second;
    return std::accumulate(l.cbegin(), l.cend(), 0, [&r](double a, auto const& b) {
        auto r_pair = r.find(b.first);
        if (r_pair == r.end()) return a;
        return a + b.second * r_pair->second;
    });
}

void QGramVector::CalculateLength() {
    length_ = std::sqrt(
        std::accumulate(q_grams_.cbegin(), q_grams_.cend(), 0, [](double a, auto const& pair) {
            return a + pair.second * pair.second;
        }));
}

}  // namespace util

#include "LevenshteinDistance.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace util {

/* Levenshtein distance computation algorithm taken from
 * https://en.wikipedia.org/wiki/Levenshtein_distance
 */
unsigned LevenshteinDistance(std::string_view l, std::string_view r) {
    size_t r_size = r.size();
    std::vector<unsigned> v0(r_size + 1);
    std::vector<unsigned> v1(r_size + 1);

    for (unsigned i = 0; i != r_size + 1; ++i) {
        v0[i] = i;
    }

    for (unsigned i = 0; i != l.size(); ++i) {
        v1[0] = i + 1;

        for (unsigned j = 0; j != r.size(); ++j) {
            assert(j + 1 < v0.size());
            unsigned del_cost = v0[j + 1] + 1;
            unsigned insert_cost = v1[j] + 1;
            unsigned substition_cost;
            if (l[i] == r[j]) {
                substition_cost = v0[j];
            } else {
                substition_cost = v0[j] + 1;
            }

            v1[j + 1] = std::min({del_cost, insert_cost, substition_cost});
        }

        std::swap(v0, v1);
    }

    return v0.back();
}

}  // namespace util

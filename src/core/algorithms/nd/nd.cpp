#include "algorithms/nd/nd.h"

#include <numeric>

#include "algorithms/nd/nd_verifier/util/vector_to_string.h"

namespace model {

using namespace algos::nd_verifier::util;

std::string ND::ToShortString() const {
    return VectorToString(GetLhsIndices()) + " -> " + VectorToString(GetRhsIndices());
}

std::string ND::ToLongString() const {
    auto cols_to_long_string = [](std::vector<Column const*> const& cols) -> std::string {
        if (cols.size() == 0) {
            return "[]";
        }

        return std::accumulate(std::next(cols.begin()), cols.end(), cols[0]->GetName(),
                               [](std::string&& a, Column const* b) {
                                   return std::move(a) + ", " + b->GetName();
                               });
    };

    return cols_to_long_string(GetLhs().GetColumns()) + " -> " +
           cols_to_long_string(GetRhs().GetColumns());
}

}  // namespace model

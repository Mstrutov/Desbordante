#include "attribute_pair.h"

namespace algos::fastod {

AttributePair::AttributePair() : AttributePair(0, 0) {}

AttributePair::AttributePair(model::ColumnIndex left, model::ColumnIndex right)
    : left(left), right(right) {}

std::string AttributePair::ToString() const {
    return "{ " + std::to_string(left + 1) + ", " + std::to_string(right + 1) + " }";
}

bool operator==(AttributePair const& x, AttributePair const& y) {
    return x.left == y.left && x.right == y.right;
}

bool operator!=(AttributePair const& x, AttributePair const& y) {
    return !(x == y);
}

bool operator<(AttributePair const& x, AttributePair const& y) {
    if (x.left != y.left) {
        return x.left < y.left;
    }

    return x.right < y.right;
}

}  // namespace algos::fastod

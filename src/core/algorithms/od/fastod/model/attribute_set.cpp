#include "attribute_set.h"

#include <sstream>

namespace algos::fastod {

std::string AttributeSet::ToString() const {
    std::stringstream result;
    result << "{";

    bool first = true;

    Iterate([&result, &first](model::ColumnIndex i) {
        if (first)
            first = false;
        else
            result << ",";

        result << i + 1;
    });

    result << "}";

    return result.str();
}

void AttributeSet::Iterate(std::function<void(model::ColumnIndex)> callback) const {
    for (model::ColumnIndex attr = FindFirst(); attr != Size(); attr = FindNext(attr)) {
        callback(attr);
    }
}

}  // namespace algos::fastod

#include "Itemset.h"

#include <algorithm>

namespace model {

void Itemset::Sort() {
    std::sort(indices_.begin(), indices_.end());
}

} // namespace model

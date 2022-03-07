#include "Itemset.h"

#include <algorithm>

void Itemset::Sort() {
    std::sort(indices_.begin(), indices_.end());
}

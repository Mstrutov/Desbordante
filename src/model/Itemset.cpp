#include "Itemset.h"

#include <algorithm>

void Itemset::sort() {
    std::sort(indices.begin(), indices.end());
}

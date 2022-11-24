#include "model/cfd.h"

// see ../algorithms/cfd/LICENSE

#include "algorithms/cfd/util/set_util.h"

[[maybe_unused]] bool IsValid(const Itemset& lhs, int rhs) {
    // Discard Variable -> Constant and Constant -> Variable CFDs
    if (rhs < 0) {
        if (lhs.size() && !Has(lhs, [](int si) -> bool { return si < 0; })) return false;
    }
    else {
        if (Has(lhs, [](int si) -> bool { return si < 0; })) return false;
    }
    return true;
}
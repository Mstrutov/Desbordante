#include "ind.h"

#include <sstream>

#include "model/table/column.h"

namespace model {

std::string IND::ToShortString() const {
    auto cc_to_short_string = [](ColumnCombination const& cc) -> std::string {
        std::stringstream ss;
        auto table_idx = cc.GetTableIndex();
        auto const& indices = cc.GetColumnIndices();
        ss << '(' << table_idx << ", [";
        for (auto it = indices.begin(); it != indices.end(); ++it) {
            if (it != indices.begin()) {
                ss << ", ";
            }
            ss << *it;
        }
        ss << "])";
        return ss.str();
    };

    std::stringstream ss;
    ss << cc_to_short_string(GetLhs()) << " -> " << cc_to_short_string(GetRhs());
    return ss.str();
}

std::string IND::ToLongString() const {
    auto cc_to_long_string = [this](ColumnCombination const& cc) -> std::string {
        std::stringstream ss;
        auto table_idx = cc.GetTableIndex();
        auto const& indices = cc.GetColumnIndices();
        ss << '(' << schemas_->at(table_idx).GetName() << ", [";
        for (auto it = indices.begin(); it != indices.end(); ++it) {
            if (it != indices.begin()) {
                ss << ", ";
            }
            ss << schemas_->at(table_idx).GetColumn(*it).GetName();
        }
        ss << "])";
        return ss.str();
    };

    std::stringstream ss;
    ss << cc_to_long_string(GetLhs()) << " -> " << cc_to_long_string(GetRhs());
    return ss.str();
}

}  // namespace model

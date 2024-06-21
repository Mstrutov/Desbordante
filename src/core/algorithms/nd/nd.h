#pragma once

#include "model/table/column.h"
#include "model/table/vertical.h"

namespace model {

using WeightType = unsigned int;

class ND {
private:
    Vertical lhs_;
    Vertical rhs_;
    WeightType weight_;

public:
    ND(Vertical lhs, Vertical rhs, WeightType weight)
        : lhs_(std::move(lhs)), rhs_(std::move(rhs)), weight_(weight) {}

    [[nodiscard]] Vertical const& GetLhs() const {
        return lhs_;
    }

    [[nodiscard]] Vertical const& GetRhs() const {
        return rhs_;
    }

    [[nodiscard]] WeightType GetWeight() const {
        return weight_;
    }

    [[nodiscard]] std::vector<ColumnIndex> GetLhsIndices() const {
        return lhs_.GetColumnIndicesAsVector();
    }

    [[nodiscard]] std::vector<ColumnIndex> GetRhsIndices() const {
        return rhs_.GetColumnIndicesAsVector();
    }

    [[nodiscard]] std::string ToShortString() const;

    [[nodiscard]] std::string ToLongString() const;
};

}  // namespace model

#pragma once
#include <string>
#include <utility>

#include "model/table/column.h"

namespace algos {
struct Correlation {
private:
    Column lhs_;
    Column rhs_;

public:
    Correlation(Column lhs, Column rhs) : lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    [[nodiscard]] Column const& GetLhs() const {
        return lhs_;
    }

    [[nodiscard]] Column::IndexType GetLhsIndex() const {
        return lhs_.GetIndex();
    }

    [[nodiscard]] std::string GetLhsName() const {
        return lhs_.GetName();
    }

    [[nodiscard]] Column const& GetRhs() const {
        return rhs_;
    }

    [[nodiscard]] Column::IndexType GetRhsInex() const {
        return rhs_.GetIndex();
    }

    [[nodiscard]] std::string GetRhsName() const {
        return rhs_.GetName();
    }

    [[nodiscard]] std::string ToString() const {
        std::string result = lhs_.GetName() + " ~ " + rhs_.GetName();
        return result;
    }
};
}  // namespace algos

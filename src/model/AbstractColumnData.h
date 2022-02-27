#pragma once

#include <string>

#include "Column.h"

namespace model {

class AbstractColumnData {
protected:
    Column const* column_;

    explicit AbstractColumnData(Column const* column) noexcept : column_(column) {}

public:
    AbstractColumnData(AbstractColumnData const& other) = default;
    AbstractColumnData& operator=(AbstractColumnData const& other) = default;
    AbstractColumnData(AbstractColumnData&& other) noexcept = default;
    AbstractColumnData& operator=(AbstractColumnData&& other) noexcept = default;

    Column const* GetColumn() const { return column_; }

    virtual std::string ToString() const = 0;
    virtual ~AbstractColumnData() = default;
};

} // namespace model


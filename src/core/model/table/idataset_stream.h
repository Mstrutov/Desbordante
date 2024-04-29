#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace model {

class IDatasetStream {
public:
    virtual std::vector<std::string> GetNextRow() = 0;
    [[nodiscard]] virtual bool HasNextRow() const = 0;
    [[nodiscard]] virtual std::size_t GetNumberOfColumns() const = 0;
    [[nodiscard]] virtual std::string GetColumnName(std::size_t index) const = 0;
    [[nodiscard]] virtual std::string GetRelationName() const = 0;
    virtual void Reset() = 0;
    virtual ~IDatasetStream() = default;

    std::vector<std::string> GetColumnNames() const {
        std::vector<std::string> column_names;
        std::size_t const number_of_columns = GetNumberOfColumns();
        column_names.reserve(number_of_columns);
        for (std::size_t i = 0; i < number_of_columns; ++i) {
            column_names.push_back(GetColumnName(i));
        }
        return column_names;
    }
};

}  // namespace model

#include "column_layout_typed_relation_data.h"

#include <easylogging++.h>

namespace model {

std::unique_ptr<ColumnLayoutTypedRelationData> ColumnLayoutTypedRelationData::CreateFrom(
        IDatasetStream& data_stream, bool is_null_eq_null) {
    auto schema = RelationalSchema::CreateFrom(data_stream);
    size_t const num_columns = data_stream.GetNumberOfColumns();

    std::vector<std::vector<std::string>> columns(num_columns);
    std::vector<std::string> row;

    /* Parsing is very similar to ColumnLayoutRelationData::CreateFrom().
     * Maybe we need column-based parsing in addition to row-based in CSVParser
     * (now IDatasetStream) */
    while (data_stream.HasNextRow()) {
        row = data_stream.GetNextRow();

        if (row.size() != num_columns) {
            LOG(WARNING) << "Unexpected number of columns for a row, skipping (expected "
                         << num_columns << ", got " << row.size() << ")";
            continue;
        }

        for (size_t index = 0; index < row.size(); ++index) {
            std::string& field = row[index];
            columns[index].push_back(std::move(field));
        }
    }

    std::vector<TypedColumnData> column_data;
    for (size_t i = 0; i < num_columns; ++i) {
        TypedColumnData typed_column_data = model::TypedColumnDataFactory::CreateFrom(
                &schema->GetColumn(i), std::move(columns[i]), is_null_eq_null);
        column_data.emplace_back(std::move(typed_column_data));
    }

    return std::make_unique<ColumnLayoutTypedRelationData>(std::move(schema),
                                                           std::move(column_data));
}

}  // namespace model

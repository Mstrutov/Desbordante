#include "ColumnLayoutTypedRelationData.h"

#include <easylogging++.h>

namespace model {

std::unique_ptr<ColumnLayoutTypedRelationData> ColumnLayoutTypedRelationData::CreateFrom(
    CSVParser& file_input, bool is_null_eq_null, int max_cols, long max_rows) {
    auto schema = std::make_unique<RelationalSchema>(file_input.GetRelationName(), is_null_eq_null);
    int num_columns = file_input.GetNumberOfColumns();

    if (max_cols > 0) {
        num_columns = std::min(num_columns, max_cols);
    }

    std::vector<std::vector<std::string>> columns(num_columns);
    int row_num = 0;
    std::vector<std::string> row;

    /* Parsing is very similar to ColumnLayoutRelationData::CreateFrom().
     * Maybe we need column-based parsing in addition to row-based in CSVParser */
    while (file_input.GetHasNext()) {
        row = file_input.ParseNext();

        if (row.empty() && num_columns == 1) {
            row.emplace_back("");
        } else if ((int)row.size() != num_columns) {
            LOG(WARNING) << "Skipping incomplete rows";
            continue;
        }

        if (max_rows <= 0 || row_num < max_rows) {
            int index = 0;
            for (std::string& field : row) {
                columns[index].push_back(std::move(field));
                index++;
                if (index >= num_columns) {
                    break;
                }
            }
        } else {
            LOG(WARNING) << "Processed " << row_num << " rows and " << max_rows - row_num
                         << " rows remain unprocessed due to `max_rows` parameter.";
            break;
        }
        row_num++;
    }

    std::vector<TypedColumnData> column_data;
    for (int i = 0; i < num_columns; ++i) {
        Column column(schema.get(), file_input.GetColumnName(i), i);
        schema->AppendColumn(std::move(column));
        TypedColumnData typed_column_data = model::TypedColumnDataFactory::CreateFrom(
            schema->GetColumn(i), std::move(columns[i]), is_null_eq_null);
        column_data.emplace_back(std::move(typed_column_data));
    }

    schema->Init();

    return std::make_unique<ColumnLayoutTypedRelationData>(std::move(schema),
                                                           std::move(column_data));
}

}  // namespace model

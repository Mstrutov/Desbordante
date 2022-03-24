//
// Created by Ilya Vologin
// https://github.com/cupertank
//
#include "ColumnLayoutRelationData.h"

#include <map>
#include <memory>
#include <utility>

#include <easylogging++.h>

std::vector<int> ColumnLayoutRelationData::GetTuple(int tuple_index) const {
    int num_columns = schema_->GetNumColumns();
    std::vector<int> tuple = std::vector<int>(num_columns);
    for (int column_index = 0; column_index < num_columns; column_index++) {
        tuple[column_index] = column_data_[column_index].GetProbingTableValue(tuple_index);
    }
    return tuple;
}

std::unique_ptr<ColumnLayoutRelationData> ColumnLayoutRelationData::CreateFrom(
    CSVParser& file_input, bool is_null_eq_null, int max_cols, long max_rows) {
    auto schema = std::make_unique<RelationalSchema>(file_input.GetRelationName(), is_null_eq_null);
    std::unordered_map<std::string, int> value_dictionary;
    int next_value_id = 1;
    const int null_value_id = -1;
    int num_columns = file_input.GetNumberOfColumns();
    if (max_cols > 0) num_columns = std::min(num_columns, max_cols);
    std::vector<std::vector<int>> column_vectors = std::vector<std::vector<int>>(num_columns);
    int row_num = 0;
    std::vector<std::string> row;

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
                if (field.empty()) {
                    column_vectors[index].push_back(null_value_id);
                } else {
                    auto location = value_dictionary.find(field);
                    int value_id;
                    if (location == value_dictionary.end()) {
                        value_dictionary[field] = next_value_id;
                        value_id = next_value_id;
                        next_value_id++;
                    } else {
                        value_id = location->second;
                    }
                    column_vectors[index].push_back(value_id);
                }
                index++;
                if (index >= num_columns) break;
            }
        } else {
            //TODO: Подумать что тут сделать
            assert(0);
        }
        row_num++;
    }

    std::vector<ColumnData> column_data;
    for (int i = 0; i < num_columns; ++i) {
        auto column = Column(schema.get(), file_input.GetColumnName(i), i);
        schema->AppendColumn(std::move(column));
        auto pli = util::PositionListIndex::CreateFor(column_vectors[i], schema->IsNullEqualNull());
        column_data.emplace_back(schema->GetColumn(i), std::move(pli));
    }

    schema->Init();

    return std::make_unique<ColumnLayoutRelationData>(std::move(schema), std::move(column_data));
}


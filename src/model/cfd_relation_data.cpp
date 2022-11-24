#include "model/cfd_relation_data.h"

// see ./LICENSE

#include <easylogging++.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <cstddef>

#include "algorithms/cfd/util/set_util.h"


size_t CFDRelationData::GetNumRows() const
{
    return data_rows_.size();
}

std::unique_ptr<CFDRelationData> CFDRelationData::CreateFrom(
    model::IDatasetStream& parser, bool is_null_eq_null,
    unsigned columns_number, unsigned tuples_number,
    double c_sample, double r_sample) {
        if (columns_number == 0 || tuples_number == 0) {
            return CFDRelationData::CreateFrom(parser, is_null_eq_null, c_sample, r_sample);
        }
        else {
            // Fields of CFDRelationData class
            auto schema = std::make_unique<RelationalSchema>(parser.GetRelationName(), is_null_eq_null);
            std::vector<Transaction> data_rows;
            boost::unordered_map<std::pair<int, std::string>, int, pairhash> item_dictionary;
            std::vector<ItemInfo> items;
            std::unordered_map<int, std::vector<int>> columns_values_dict;
            int unique_elems_number = 1;
        
            unsigned num_columns = parser.GetNumberOfColumns();
            std::vector<std::string> line;
            num_columns = std::min(num_columns, columns_number);
            std::vector<std::string> string_row(num_columns);
            while (parser.HasNextRow() && data_rows.size() < tuples_number) {
                line = parser.GetNextRow();
                std::vector<int> int_row(num_columns);
                for (size_t k = 0; k < num_columns; k++) {
                    string_row[k] = line[k];
                }
                int it;
                for (size_t i = 0; i < num_columns; i++) {
                    auto ptr = item_dictionary.find(std::make_pair(i, string_row[i]));
                    if (ptr != item_dictionary.end()) {
                        it = ptr->second;
                    }
                    else {
                        items.emplace_back(string_row[i], i);
                        columns_values_dict[(size_t)i].push_back(unique_elems_number);
                        item_dictionary[std::make_pair(i, string_row[i])] = unique_elems_number;
                        it = unique_elems_number++;
                    }
                    int_row[i] = it;
                    items[it - 1].frequency++;
                }
                data_rows.push_back(int_row);
            }

            std::vector<CFDColumnData> column_data;
            for (size_t i = 0; i < num_columns; ++i) {
                auto column = Column(schema.get(), parser.GetColumnName(i), i);
                schema->AppendColumn(std::move(column));
                column_data.emplace_back(schema->GetColumn(i), columns_values_dict[i]);
            }
            schema->Init();
            return std::make_unique<CFDRelationData>(std::move(schema), std::move(column_data), 
            std::move(data_rows), std::move(item_dictionary), std::move(items));
        }
    }

std::unique_ptr<CFDRelationData> CFDRelationData::CreateFrom(
    model::IDatasetStream& file_input, bool is_null_eq_null, double c_sample, double rSample) {
    // Fields of CFDRelationData class
    auto schema = std::make_unique<RelationalSchema>(file_input.GetRelationName(), is_null_eq_null);
    std::vector<Transaction> data_rows;
    boost::unordered_map<std::pair<int, std::string>, int, pairhash> item_dictionary;
    std::vector<ItemInfo> items;
    std::unordered_map<int, std::vector<int>> columns_values_dict;
    int unique_elems_number = 1;
    std::random_device rd; // only used once to initialise (seed) engine
    std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<double> uni(0.0,1.0); // guaranteed unbiased
    int num_columns = (int)file_input.GetNumberOfColumns();
	std::vector<std::string> line;
    std::vector<int> columns_numbers_list = Range(0, num_columns);
    int size = columns_numbers_list.size() * c_sample;
    Shuffle(columns_numbers_list);
    columns_numbers_list = std::vector<int>(columns_numbers_list.begin(),
                                            columns_numbers_list.begin() + size);
    std::sort(columns_numbers_list.begin(), columns_numbers_list.end());
    while (file_input.HasNextRow()) {
        if (uni(rng) < rSample) {
        line = file_input.GetNextRow();
        std::vector<int> row(size);
        int j = 0;
        int it;
        for (size_t i = 0; i < line.size(); i++) {
            if (std::binary_search(columns_numbers_list.begin(), columns_numbers_list.end(), i)) {
                auto ptr = item_dictionary.find(std::make_pair(j, line[i]));
                if (ptr != item_dictionary.end()) {
                    it = ptr->second;
                }
                else {
                    items.emplace_back(line[i], j);
                    columns_values_dict[j].push_back(unique_elems_number);
                    item_dictionary[std::make_pair(j, line[i])] = unique_elems_number;
                    it = unique_elems_number++;
                }
                items[it - 1].frequency++;
                row[j] = it;
                j++;
            }
        }
        if (j > 0) {
            data_rows.push_back(row);
        }
    }
    }

    std::vector<CFDColumnData> column_data;
    for (int i = 0; i < num_columns; ++i) {
        auto column = Column(schema.get(), file_input.GetColumnName(i), i);
        schema->AppendColumn(std::move(column));
        column_data.emplace_back(schema->GetColumn(i), columns_values_dict[i]);
    }
    schema->Init();
    return std::make_unique<CFDRelationData>(std::move(schema), std::move(column_data), 
    std::move(data_rows), std::move(item_dictionary), std::move(items));
}

unsigned CFDRelationData::Size() const {
    return data_rows_.size();
}

unsigned CFDRelationData::GetAttrsNumber() const {
    return column_data_.size();
}

unsigned CFDRelationData::GetItemsNumber() const {
    return items_.size();
}

const Transaction& CFDRelationData::GetRow(int row) const {
    return data_rows_.at(row);
}

void CFDRelationData::SetRow(int row_index, const Transaction& row) {
    for (size_t i = 0; i < row.size(); i++) {
        if (data_rows_[row_index][i] != row[i]) {
            items_[data_rows_[row_index][i] - 1].frequency--;
            items_[row[i] - 1].frequency++;
        }
    }
    data_rows_[row_index] = row;
}

const std::vector<int>& CFDRelationData::GetDomainOfItem(int item) const {
    return column_data_.at(items_[item - 1].attribute).GetValues();
}

const std::vector<int>& CFDRelationData::GetDomain(int attr) const {
    return column_data_.at(attr).GetValues();
}

std::string CFDRelationData::GetAttrName(int index) const
{
    return GetSchema()->GetColumn(index)->GetName();
}

int CFDRelationData::GetAttr(const std::string& s) const {
    for (size_t i = 0; i < GetNumColumns(); i++) {
        const std::string ai = GetAttrName(i);
        if (ai == s) {
            return i;
        }
    }
    return -1;
}

int CFDRelationData::GetItem(int attr, const std::string& str_value) const {
    return item_dictionary_.at(std::make_pair(attr, str_value));
}

void CFDRelationData::Sort() {
    std::sort(data_rows_.begin(), data_rows_.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return std::lexicographical_compare(a.begin(), b.begin(), a.end(), b.end());
        });
}

void CFDRelationData::ToFront(const SimpleTidList& tids) {
    for (size_t i = 0; i < tids.size(); i++) {
        std::swap(data_rows_[i], data_rows_[tids[i]]);
    }
}

int CFDRelationData::Frequency(int i) const {
    return items_[i - 1].frequency;
}

int CFDRelationData::GetAttrIndex(int item_index) const {
    if (item_index > 0) return items_[item_index - 1].attribute;
    return -1 - item_index;
}

std::string CFDRelationData::GetStringFormat(char delim) const {
    std::string file;
    for (size_t ai = 0; ai < GetNumColumns(); ai++) {
        const std::string attr = GetAttrName((int)ai);
        file += attr;
        if (ai < GetNumColumns() - 1)
            file += delim;
        else
            file += '\n';
    }
    for (const auto& row : data_rows_) {
        for (size_t ri = 0; ri < row.size(); ri++) {
            const auto& item = row[ri];
            file += items_[item-1].value;
            if (ri < row.size()-1)
                file += delim;
            else
                file += '\n';
        }
    }
    return file;
}

std::string CFDRelationData::GetStringFormat(const SimpleTidList& subset, char delim) const {
    std::string result;
    for (size_t ai = 0; ai < GetNumColumns(); ai++) {
        const auto& attr = GetAttrName((int)ai);
        result += attr;
        if (ai < GetNumColumns() - 1)
            result += delim;
        else
            result += '\n';
    }
    for (int i : subset) {
        auto& row = data_rows_[i];
        for (size_t ri = 0; ri < row.size(); ri++) {
            const auto& item = row[ri];
            result += items_[item - 1].value;
            if (ri < row.size() - 1)
                result +=  delim;
            else
                result += '\n';
        }
    }
    return result;
}


const std::string& CFDRelationData::GetValue(int i) const {
    return items_[i - 1].value;
}

// Receives elements from database and returns its attributes sorted vector
std::vector<int> CFDRelationData::GetAttrVector(const Itemset& items) const {
    std::vector<int> attrs;
    attrs.reserve(items.size());
    for (int i : items) {
        if (i > 0)
            attrs.push_back(GetAttrIndex(i));
        else
            attrs.push_back(-1 - i);
    }
    std::sort(attrs.begin(), attrs.end());
    return attrs;
}

std::vector<int> CFDRelationData::GetAttrVectorItems(const Itemset& items) const {
    std::vector<int> attrs;
    attrs.reserve(items.size());
    for (int i : items) {
        if (i > 0)
            attrs.push_back(-1 - GetAttrIndex(i));
        else
            attrs.push_back(i);
    }
    std::sort(attrs.begin(), attrs.end());
    return attrs;
}

#include "hashed_column_store.h"

#include "model/table/column.h"

namespace algos::faida {

std::filesystem::path HashedColumnStore::PrepareDirNext(std::filesystem::path dir,
                                                        TableIndex table_idx) {
    for (Column const& column : schema_->GetColumns()) {
        std::string file_name;
        file_name += std::to_string(table_idx);
        file_name += "_";
        file_name += std::to_string(column.GetIndex());
        file_name += ".bin";

        std::filesystem::path column_file = dir / file_name;
        column_files_[column.GetIndex()] = std::move(column_file);
    }
    return dir;
}

void HashedColumnStore::WriteColumnsAndSample(model::IDatasetStream& data_stream) {
    bool is_writing_any_column = false;
    std::vector<std::ofstream> column_files_out(schema_->GetNumColumns());

    ColumnIndex column_idx = 0;
    for (std::filesystem::path const& column_file : column_files_) {
        column_files_out[column_idx++].open(column_file, std::ios::binary);
        is_writing_any_column = true;
    }

    std::vector<std::vector<std::string>> rows_to_sample;
    std::vector<std::optional<size_t>> constant_col_hashes(schema_->GetNumColumns(), std::nullopt);
    std::vector sampled_col_values(
            schema_->GetNumColumns(),
            std::unordered_set<size_t>(sample_goal_ < 0 ? 1000 : sample_goal_));

    int const bufsize = read_buff_size_;
    std::vector buf(schema_->GetNumColumns(), std::vector<size_t>(bufsize));

    int row_counter = 0;
    while (data_stream.HasNextRow()) {
        std::vector<std::string> row = data_stream.GetNextRow();
        if (row.empty() || row.size() != schema_->GetNumColumns()) {
            continue;
        }

        bool is_sample_complete = true;
        bool row_has_unseen_value = false;
        for (ColumnIndex col_idx = 0; col_idx < schema_->GetNumColumns(); col_idx++) {
            std::string const& value = row.at(col_idx);
            size_t value_hash = this->Hash(value);

            buf[col_idx][row_counter % bufsize] = value_hash;
            if (row_counter % bufsize == bufsize - 1) {
                column_files_out[col_idx].write(reinterpret_cast<char*>(buf[col_idx].data()),
                                                bufsize * sizeof(size_t));
            }

            if (row_counter == 0) {
                // Assume all columns are constant initially
                constant_col_hashes[col_idx] = value_hash;
            } else if (constant_col_hashes[col_idx].has_value() &&
                       constant_col_hashes[col_idx].value() != value_hash) {
                constant_col_hashes[col_idx].reset();
            }

            if (value_hash != null_hash_) {
                std::unordered_set<size_t>& sampled_vals = sampled_col_values[col_idx];
                bool const should_sample =
                        sample_goal_ < 0 ||
                        sampled_vals.size() < static_cast<unsigned>(sample_goal_);

                is_sample_complete &= !should_sample;
                if (should_sample && sampled_vals.insert(value_hash).second) {
                    row_has_unseen_value = true;
                }
            }
        }

        if (row_has_unseen_value) {
            rows_to_sample.emplace_back(std::move(row));
        }

        row_counter++;

        if (!is_writing_any_column && is_sample_complete) {
            break;
        }
    }

    if (row_counter == 0) {
        throw std::runtime_error("Got an empty file: IND mining is meaningless.");
    }

    // Write data from the last buffer
    for (ColumnIndex col_idx = 0; col_idx < schema_->GetNumColumns(); col_idx++) {
        column_files_out[col_idx].write(reinterpret_cast<char*>(buf[col_idx].data()),
                                        (row_counter % bufsize) * sizeof(size_t));
    }

    WriteSample(rows_to_sample);

    for (ColumnIndex col_idx = 0; col_idx < constant_col_hashes.size(); col_idx++) {
        if (constant_col_hashes[col_idx].has_value()) {
            size_t const col_hash = constant_col_hashes[col_idx].value();
            if (col_hash == null_hash_) {
                column_properties_[col_idx] = ColumnProperty::kNullConstant;
            } else {
                column_properties_[col_idx] = ColumnProperty::kConstant;
            }
        }
    }

    for (std::ofstream& col_file : column_files_out) {
        col_file.close();
    }
}

std::unique_ptr<IRowIterator> HashedColumnStore::GetRows(
        std::unordered_set<ColumnIndex> const& columns) const {
    std::vector<std::optional<std::ifstream>> hashed_col_streams(schema_->GetNumColumns());

    for (ColumnIndex col_idx : columns) {
        hashed_col_streams[col_idx] = std::ifstream(column_files_[col_idx], std::ios::binary);
    }

    return std::make_unique<RowIterator>(std::move(hashed_col_streams));
}

std::unique_ptr<AbstractColumnStore> HashedColumnStore::CreateFrom(
        std::string const& dataset_name, TableIndex table_idx, model::IDatasetStream& input_data,
        int sample_goal, size_t null_hash) {
    auto store = std::make_unique<HashedColumnStore>(
            HashedColumnStore(input_data.GetNumberOfColumns(), sample_goal, null_hash));
    store->LoadData(dataset_name, table_idx, input_data);

    return store;
}

HashedColumnStore::RowIterator::~RowIterator() {
    for (std::optional<std::ifstream>& column_stream : hashed_col_streams_) {
        if (column_stream.has_value()) {
            column_stream->close();
        }
    }
}

bool HashedColumnStore::RowIterator::HasNextBlock() {
    if (!has_next_) return false;

    curr_block_size_ = default_block_size_;
    Block row_hashes_inv(hashed_col_streams_.size());

    ColumnIndex col_idx = 0;
    for (std::optional<std::ifstream>& column_stream : hashed_col_streams_) {
        if (column_stream.has_value()) {
            row_hashes_inv[col_idx] = AlignedVector(default_block_size_);
            if (!column_stream->read(reinterpret_cast<char*>(row_hashes_inv[col_idx]->data()),
                                     default_block_size_ * sizeof(size_t))) {
                has_next_ = false;
                curr_block_size_ = column_stream->gcount() / sizeof(size_t);
            }
        }
        ++col_idx;
    }

    curr_block_ = std::move(row_hashes_inv);
    return true;
}

IRowIterator::Block const& HashedColumnStore::RowIterator::GetNextBlock() {
    return curr_block_;
}

}  // namespace algos::faida

#include "abstract_column_store.h"

#include <fstream>

#include "algorithms/ind/faida/hashing/hashing.h"
#include "model/table/column.h"

namespace algos::faida {

void AbstractColumnStore::LoadData(std::string const& dataset_name, TableIndex table_idx,
                                   model::IDatasetStream& input_data) {
    size_t const num_columns = input_data.GetNumberOfColumns();
    if (num_columns == 0) {
        throw std::runtime_error("Got an empty file: IND mining is meaningless.");
    }

    schema_ = RelationalSchema::CreateFrom(input_data);
    column_properties_ =
            std::vector<ColumnProperty>(input_data.GetNumberOfColumns(), ColumnProperty::kOrdinary);

    std::filesystem::path dir = PrepareDir(dataset_name, table_idx);

    WriteColumnsAndSample(input_data);
}

std::filesystem::path AbstractColumnStore::PrepareDir(std::string const& dataset_name,
                                                      TableIndex table_idx) {
    namespace fs = std::filesystem;
    fs::path sample_location_dir =
            fs::current_path() / dir_name_ / dataset_name / schema_->GetName();

    fs::create_directories(sample_location_dir);

    std::string sample_file_name;
    sample_file_name += std::to_string(table_idx);
    sample_file_name += "-sample.bin";
    sample_file_ = sample_location_dir / sample_file_name;

    return PrepareDirNext(sample_location_dir, table_idx);
}

void AbstractColumnStore::WriteSample(std::vector<std::vector<std::string>> const& rows) {
    sample_size_ = rows.size();
    std::ofstream sample_stream(sample_file_, std::ios::binary);

    std::vector<size_t> row_hashes(schema_->GetNumColumns());
    for (std::vector<std::string> const& row : rows) {
        ColumnIndex col_idx = 0;
        for (std::string const& value : row) {
            size_t const value_hash = Hash(value);
            row_hashes[col_idx++] = value_hash;
        }
        sample_stream.write(reinterpret_cast<char*>(row_hashes.data()),
                            sizeof(size_t) * row_hashes.size());
    }

    sample_stream.close();
}

AbstractColumnStore::HashedTableSample AbstractColumnStore::ReadSample() const {
    std::ifstream sample_stream(sample_file_, std::ios::binary);

    HashedTableSample hashed_sample(sample_size_);
    for (unsigned row_idx = 0; row_idx < sample_size_; row_idx++) {
        std::vector<size_t> row_hashes(schema_->GetNumColumns());
        sample_stream.read(reinterpret_cast<char*>(row_hashes.data()),
                           sizeof(size_t) * row_hashes.size());

        hashed_sample[row_idx] = std::move(row_hashes);
    }

    sample_stream.close();
    return hashed_sample;
}

}  // namespace algos::faida

#include "algorithms/md/hymd/indexes/dictionary_compressor.h"

#include <cstddef>

#include <easylogging++.h>

#include "model/index.h"

namespace algos::hymd::indexes {

DictionaryCompressor::DictionaryCompressor(std::size_t attribute_num) : plis_(attribute_num) {}

void DictionaryCompressor::AddRecord(std::vector<std::string> record) {
    std::size_t const record_size = record.size();
    std::size_t const expected_size = plis_.size();
    if (record_size != expected_size) {
        LOG(WARNING) << "Unexpected number of columns for a record, skipping (expected "
                     << expected_size << ", got " << record_size
                     << "). Records processed so far: " << records_processed_ << ".";
        return;
    }
    ++records_processed_;
    CompressedRecord rec;
    rec.reserve(expected_size);
    for (model::Index i = 0; i < expected_size; ++i) {
        rec.push_back(plis_[i].AddNextValue(record[i]));
    }
    records_.push_back(std::move(rec));
}

std::unique_ptr<DictionaryCompressor> DictionaryCompressor::CreateFrom(
        model::IDatasetStream& stream /*, indices*/) {
    std::size_t const columns = stream.GetNumberOfColumns();
    auto compressor = std::make_unique<DictionaryCompressor>(columns);
    if (!stream.HasNextRow())
        throw std::runtime_error("MD mining is meaningless on empty dataset!");
    while (stream.HasNextRow()) compressor->AddRecord(stream.GetNextRow());
    return compressor;
}

}  // namespace algos::hymd::indexes

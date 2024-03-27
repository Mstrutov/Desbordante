#pragma once

#include "algorithms/md/hymd/indexes/dictionary_compressor.h"

namespace algos::hymd::indexes {

class RecordsInfo {
    std::shared_ptr<DictionaryCompressor const> const left_compressor_;
    std::shared_ptr<DictionaryCompressor const> const right_compressor_;

public:
    [[nodiscard]] bool OneTableGiven() const noexcept {
        return left_compressor_ == right_compressor_;
    }

    [[nodiscard]] DictionaryCompressor const& GetLeftCompressor() const noexcept {
        return *left_compressor_;
    }

    [[nodiscard]] DictionaryCompressor const& GetRightCompressor() const noexcept {
        return *right_compressor_;
    }

    RecordsInfo(std::shared_ptr<DictionaryCompressor> left_compressor,
                std::shared_ptr<DictionaryCompressor> right_compressor) noexcept
        : left_compressor_(std::move(left_compressor)),
          right_compressor_(std::move(right_compressor)) {}

    static std::unique_ptr<RecordsInfo> CreateFrom(model::IDatasetStream& left_table) {
        std::shared_ptr<DictionaryCompressor> left_compressor =
                DictionaryCompressor::CreateFrom(left_table);
        return std::make_unique<RecordsInfo>(left_compressor, left_compressor);
    }

    static std::unique_ptr<RecordsInfo> CreateFrom(model::IDatasetStream& left_table,
                                                   model::IDatasetStream& right_table) {
        return std::make_unique<RecordsInfo>(DictionaryCompressor::CreateFrom(left_table),
                                             DictionaryCompressor::CreateFrom(right_table));
    }
};

}  // namespace algos::hymd::indexes

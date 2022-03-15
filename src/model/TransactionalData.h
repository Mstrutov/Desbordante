#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "CSVParser.h"
#include "Itemset.h"
#include "TransactionalInputFormat.h"

namespace model {

class TransactionalData {
private:
    std::vector<std::string> item_universe_;
    std::unordered_map<unsigned, Itemset> transactions_;

    static std::unique_ptr<TransactionalData> CreateFromSingular(CSVParser& file_input,
                                                                 unsigned tid_col_index = 0,
                                                                 unsigned item_col_index = 1);

    static std::unique_ptr<TransactionalData> CreateFromTabular(CSVParser& file_input,
                                                                bool has_tid);

    TransactionalData(std::vector<std::string> item_universe,
                      std::unordered_map<unsigned, Itemset> transactions)
        : item_universe_(std::move(item_universe)), transactions_(std::move(transactions)) {}

public:
    TransactionalData() = delete;

    TransactionalData(TransactionalData const&) = delete;
    TransactionalData& operator=(TransactionalData const&) = delete;

    TransactionalData(TransactionalData&& other) = default;
    TransactionalData& operator=(TransactionalData&& other) = default;

    std::vector<std::string> const& GetItemUniverse() const noexcept { return item_universe_; }
    std::unordered_map<unsigned, Itemset> const& GetTransactions() const noexcept {
        return transactions_;
    }

    size_t GetUniverseSize() const noexcept { return item_universe_.size(); }
    size_t GetNumTransactions() const noexcept { return transactions_.size(); }

    static std::unique_ptr<TransactionalData> CreateFrom(CSVParser& file_input,
                                                         InputFormat const& input_type);
};

}  // namespace model

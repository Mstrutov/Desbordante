#include "TransactionalData.h"

#include <cassert>
#include <unordered_map>

namespace model {

std::unique_ptr<TransactionalData> TransactionalData::CreateFrom(CSVParser& file_input,
                                                                 InputFormat const& input_type) {
    if (typeid(input_type) == typeid(Singular)) {
        return TransactionalData::CreateFromSingular(file_input, input_type.tid_column_index(),
                                                     input_type.item_column_index());
    } else if (typeid(input_type) == typeid(Tabular)) {
        return TransactionalData::CreateFromTabular(file_input, input_type.tid_presence());
    } else {
        throw std::logic_error("This input type is not maintained yet");
    }
}

std::unique_ptr<TransactionalData> TransactionalData::CreateFromSingular(CSVParser& file_input,
                                                                         unsigned tid_col_index,
                                                                         unsigned item_col_index) {
    std::vector<std::string> item_universe;
    std::unordered_map<std::string, unsigned> item_universe_set;
    std::unordered_map<unsigned, Itemset> transactions;
    unsigned latest_item_id = 0;

    assert(file_input.GetNumberOfColumns() >
           static_cast<int>(std::max(tid_col_index, item_col_index)));

    while (file_input.GetHasNext()) {
        std::vector<std::string> const row = file_input.ParseNext();
        if (row.empty()) {
            continue;
        }

        unsigned const tid = std::stoi(row[tid_col_index]);
        std::string const& item_name = row[item_col_index];
        unsigned item_id = latest_item_id;

        auto const [item_iter, was_inserted] = item_universe_set.insert({item_name, item_id});
        if (was_inserted) {
            // TODO(alexandrsmirn) попробовать избежать этого добавления.
            item_universe.push_back(item_name);
            ++latest_item_id;
        } else {
            // if this item already exists in the universe, set the old item id
            item_id = item_iter->second;
        }

        auto transaction = transactions.find(tid);
        if (transaction == transactions.end()) {
            Itemset items;
            items.AddItemId(item_id);
            transactions.insert({tid, std::move(items)});
        } else {
            transaction->second.AddItemId(item_id);
        }
    }

    // sort items in each transaction
    for (auto& [tid, items] : transactions) {
        items.Sort();
    }

    return std::make_unique<TransactionalData>(std::move(item_universe), std::move(transactions));
}

std::unique_ptr<TransactionalData> TransactionalData::CreateFromTabular(CSVParser& file_input,
                                                                        bool has_tid) {
    std::vector<std::string> item_universe;
    std::unordered_map<std::string, unsigned> item_universe_set;
    std::unordered_map<unsigned, Itemset> transactions;
    unsigned latest_item_id = 0;
    unsigned tid = 0;

    while (file_input.GetHasNext()) {
        std::vector<std::string> const row = file_input.ParseNext();
        if (row.empty()) {
            continue;
        }

        auto row_iter = row.begin();
        if (has_tid) {
            tid = std::stoi(*row_iter);
            row_iter++;
        }

        Itemset items;
        for (; row_iter != row.end(); ++row_iter) {
            std::string const& item_name = *row_iter;
            if (item_name.empty()) {
                continue;
            }
            unsigned item_id = latest_item_id;

            auto const item_insertion_result = item_universe_set.insert({item_name, item_id});
            if (item_insertion_result.second) {
                // TODO(alexandrsmirn) попробовать избежать этого добавления.
                item_universe.push_back(item_name);
                ++latest_item_id;
            } else {
                // if this item already exists in the universe, set the old item id
                item_id = item_insertion_result.first->second;
            }
            items.AddItemId(item_id);
        }

        items.Sort();
        transactions.insert({tid, std::move(items)});
        if (!has_tid) {
            ++tid;
        }
    }

    return std::make_unique<TransactionalData>(std::move(item_universe), std::move(transactions));
}

} // namespace model

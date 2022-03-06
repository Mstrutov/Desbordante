#include "TransactionalData.h"

#include <cassert>
#include <unordered_map>

std::unique_ptr<TransactionalData>
TransactionalData::createFrom(CSVParser& fileInput, InputFormat const& inputFormat) {
    if (typeid(inputFormat) == typeid(Singular)) {
        return TransactionalData::createFromTwoColumns(fileInput, inputFormat.tid_column_index(),
                                                                  inputFormat.item_column_index());
    } else if (typeid(inputFormat) == typeid(Tabular)) {
        return TransactionalData::createFromItemsetRows(fileInput, inputFormat.tid_presence());
    } else {
        throw std::logic_error("This input type is not maintained yet");
    }
}

std::unique_ptr<TransactionalData>
TransactionalData::createFromTwoColumns(CSVParser& fileInput, unsigned tidColumn,
                                                              unsigned itemColumn) {
    std::vector<std::string> itemUniverse;
    std::unordered_map<std::string, unsigned> itemUniverseSet;
    std::unordered_map<unsigned, Itemset> transactions;
    unsigned lastItemID = 0;

    assert(fileInput.GetNumberOfColumns() > static_cast<int>(std::max(tidColumn, itemColumn)));
    while (fileInput.GetHasNext()) {
        std::vector<std::string> const row = fileInput.ParseNext();
        if (row.empty()) {
            continue;
        }

        unsigned const transactionId = std::stoi(row[tidColumn]);
        std::string const& itemName = row[itemColumn];
        unsigned itemID = lastItemID;

        auto const [itemIter, wasInserted] = itemUniverseSet.insert({itemName, itemID});
        if (wasInserted) {
            itemUniverse.push_back(itemName);   //TODO лишняя трата памяти
            ++lastItemID;
        } else {
            //if there is such item in the universe, get its itemID
            itemID = itemIter->second;
        }

        auto transaction = transactions.find(transactionId);
        if (transaction == transactions.end()) {
            Itemset items;
            items.addItemID(itemID);
            transactions.insert({transactionId, std::move(items)});
        } else {
            transaction->second.addItemID(itemID);
        }
    }

    //sort items in each transaction
    for (auto & [tID, items] : transactions) {
        items.sort();
    }

    return std::make_unique<TransactionalData>(std::move(itemUniverse), std::move(transactions));
}

std::unique_ptr<TransactionalData>
TransactionalData::createFromItemsetRows(CSVParser &fileInput, bool hasTransactionID) {
    std::vector<std::string> itemUniverse;
    std::unordered_map<std::string, unsigned> itemUniverseSet;
    std::unordered_map<unsigned, Itemset> transactions;
    unsigned lastItemID = 0;
    unsigned transactionId = 0;

    while (fileInput.GetHasNext()) {
        std::vector<std::string> const row = fileInput.ParseNext();
        if (row.empty()) {
            continue;
        }

        auto rowIterator = row.begin();
        if (hasTransactionID)  {
            transactionId = std::stoi(*rowIterator);
            rowIterator++;
        }

        Itemset items;
        for (; rowIterator != row.end(); ++rowIterator) {
            std::string const& itemName = *rowIterator;
            if (itemName.empty()) {
                continue;
            }
            unsigned itemID = lastItemID;

            auto const itemInsertionResult = itemUniverseSet.insert({itemName, itemID});
            if (itemInsertionResult.second) {
                itemUniverse.push_back(itemName);   //TODO лишняя трата памяти
                ++lastItemID;
            } else {                                        //if there is item in the universe
                itemID = itemInsertionResult.first->second; //set old item id
            }
            items.addItemID(itemID);
        }

        items.sort();
        transactions.insert({transactionId, std::move(items)});
        if (!hasTransactionID) {
            ++transactionId;
        }
    }

    return std::make_unique<TransactionalData>(std::move(itemUniverse), std::move(transactions));
}

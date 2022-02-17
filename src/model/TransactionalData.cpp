#include "TransactionalData.h"

#include <cassert>
#include <unordered_map>

std::unique_ptr<TransactionalData> TransactionalData::createFrom(CSVParser &fileInput,
                                                                 TransactionalInputFormat inputFormat,
                                                                 bool hasTransactionID) {
    switch (inputFormat) {
        case TransactionalInputFormat::TwoColumns:
            return createFromTwoColumns(fileInput);
        case TransactionalInputFormat::ItemsetRows:
            return createFromItemsetRows(fileInput, hasTransactionID);
        default:
            return createFromTwoColumns(fileInput);
    }
}

std::unique_ptr<TransactionalData> TransactionalData::createFromTwoColumns(CSVParser& fileInput) {
    std::vector<std::string> itemUniverse;
    std::unordered_map<std::string, unsigned> itemUniverseSet;
    std::unordered_map<unsigned, Itemset> transactions;
    unsigned lastItemID = 0;

    assert(fileInput.getNumberOfColumns() >= 2);
    while (fileInput.getHasNext()) {
        std::vector<std::string> const row = fileInput.parseNext();
        //TODO what if there is an incomplete row?

        if (row.empty()) {
            continue;
        }

        unsigned const transactionId = std::stoi(row[0]);
        std::string const& itemName = row[1];
        unsigned itemID = lastItemID;

        auto const itemInsertionResult = itemUniverseSet.insert({itemName, itemID});
        if (itemInsertionResult.second) {
            itemUniverse.push_back(itemName);   //TODO лишняя трата памяти
            ++lastItemID;
        } else {                                        //if there is item in the universe
            itemID = itemInsertionResult.first->second; //set old item id
        }

        auto transaction = transactions.find(transactionId);
        if (transaction == transactions.end()) {
            Itemset items;
            items.addItemID(itemID);
            transactions.insert({transactionId, std::move(items)});
        } else {
            transactions[transactionId].addItemID(itemID);
        }
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

    while (fileInput.getHasNext()) {
        std::vector<std::string> const row = fileInput.parseNext();

        if (row.empty()) { //TODO empty transactions?
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
                break;
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
        transactions.insert({transactionId, std::move(items)});
        if (!hasTransactionID) {
            ++transactionId;
        }
    }

    return std::make_unique<TransactionalData>(std::move(itemUniverse), std::move(transactions));
}

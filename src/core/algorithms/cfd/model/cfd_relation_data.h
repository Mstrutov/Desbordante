#pragma once

#include <cstddef>
#include <utility>

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

#include "cfd_column_data.h"
#include "cfd_types.h"
#include "model/table/idataset_stream.h"
#include "model/table/relation_data.h"

// see algorithms/cfd/LICENSE

namespace algos::cfd {

// Data presentation class that CFDDiscovery uses.
class CFDRelationData : public AbstractRelationData<CFDColumnData> {
private:
    using ItemDictionary = boost::unordered_map<std::pair<int, std::string>, int, PairHash>;
    using ColumnesValuesDict = std::unordered_map<AttributeIndex, std::vector<int>>;

    // ItemInfo contains info about one elem in the table.
    struct ItemInfo {
        ItemInfo() = default;

        ItemInfo(std::string val, AttributeIndex attr)
            : value(std::move(val)), attribute(attr), frequency(0) {}

        std::string value;
        AttributeIndex attribute;
        unsigned frequency;
    };

    // array of data represented as rows of integers
    std::vector<Transaction> data_rows_;
    // maps a value string for the given attribute index to corresponding item id
    boost::unordered_map<std::pair<int, std::string>, int, PairHash> item_dictionary_;
    std::vector<ItemInfo> items_;

    static void AddNewItemsInFullTable(ItemDictionary &, ColumnesValuesDict &,
                                       std::vector<ItemInfo> &, std::vector<std::string> const &,
                                       std::vector<int> &, std::vector<Transaction> &, int &,
                                       size_t);

    static void AddNewItemsInPartialTable(ItemDictionary &, ColumnesValuesDict &,
                                          std::vector<ItemInfo> &, std::vector<std::string> const &,
                                          std::vector<int> const &, std::vector<Transaction> &,
                                          int &, int);

public:
    unsigned Size() const;
    unsigned GetAttrsNumber() const;
    unsigned GetItemsNumber() const;
    size_t GetNumRows() const override;
    Transaction const &GetRow(unsigned) const;
    std::string GetStringFormat(char delim = ' ') const;
    std::string GetStringFormat(SimpleTIdList const &subset, char delim = ' ') const;
    void Sort();
    void ToFront(SimpleTIdList const &);
    void SetRow(int row_index, Transaction const &row);
    AttributeIndex GetAttrIndex(int item_index) const;
    unsigned Frequency(int i) const;
    std::string const &GetValue(int i) const;
    std::vector<int> const &GetDomainOfItem(int) const;
    std::vector<int> const &GetDomain(unsigned attr) const;
    std::vector<int> GetAttrVector(Itemset const &) const;
    std::vector<int> GetAttrVectorItems(Itemset const &) const;
    std::string GetAttrName(int index) const;
    int GetAttr(std::string const &) const;
    int GetItem(int, std::string const &) const;

    static std::unique_ptr<CFDRelationData> CreateFrom(model::IDatasetStream &file_input,
                                                       double c_sample = 1, double r_sample = 1);
    static std::unique_ptr<CFDRelationData> CreateFrom(model::IDatasetStream &file_input,
                                                       size_t columns_number, size_t tuples_number,
                                                       double c_sample = 1, double r_sample = 1);

    CFDRelationData(std::unique_ptr<RelationalSchema> schema,
                    std::vector<CFDColumnData> column_data, std::vector<Transaction> data,
                    boost::unordered_map<std::pair<int, std::string>, int, PairHash> item_dict,
                    std::vector<ItemInfo> items)
        : AbstractRelationData(std::move(schema), std::move(column_data)),
          data_rows_(std::move(data)),
          item_dictionary_(std::move(item_dict)),
          items_(std::move(items)) {}
};
}  // namespace algos::cfd

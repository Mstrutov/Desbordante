#pragma once

// see ./LICENSE

#include <utility>
#include <cstddef>

#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

#include "model/relation_data.h"
#include "model/cfd_column_data.h"
#include "model/cfd_types.h"
#include "model/idataset_stream.h"

// Data presentation class that CFDDiscovery uses.
class CFDRelationData : public AbstractRelationData<CFDColumnData>
{
private:
    struct ItemInfo {
        ItemInfo() = default;
        ItemInfo(std::string  val, int attr) : value(std::move(val)), attribute(attr), frequency(0) { }
        std::string value;
        int attribute;
        int frequency;
    };

    // array of data represented as rows of integers
    std::vector<Transaction> data_rows_;
    // maps a value string for the given attribute index to corresponding item id
    boost::unordered_map<std::pair<int, std::string>, int, pairhash> item_dictionary_;
    std::vector<ItemInfo> items_;

public:
	unsigned Size() const;
	unsigned GetAttrsNumber() const;
    unsigned GetItemsNumber() const;
	size_t GetNumRows() const override;
    const Transaction& GetRow(int) const;
    std::string GetStringFormat(char delim= ' ') const;
    std::string GetStringFormat(const SimpleTidList& subset, char delim= ' ') const;
    void Sort();
    void ToFront(const SimpleTidList&);
    void SetRow(int row_index, const Transaction& row);
	int GetAttrIndex(int item_index) const;
	int Frequency(int i) const;
	const std::string& GetValue(int i) const;
    const std::vector<int>& GetDomainOfItem(int) const;
	const std::vector<int>& GetDomain(int attr) const;
	std::vector<int> GetAttrVector(const Itemset&) const;
    std::vector<int> GetAttrVectorItems(const Itemset&) const;
    std::string GetAttrName(int index) const;
    int GetAttr(const std::string&) const;
    int GetItem(int, const std::string&) const;

    static std::unique_ptr<CFDRelationData> CreateFrom(model::IDatasetStream& file_input,
                                                       bool is_null_eq_null,
                                                       double c_sample = 1,
                                                       double rSample = 1);
    static std::unique_ptr<CFDRelationData> CreateFrom(model::IDatasetStream& file_input,
                                                       bool is_null_eq_null,
                                                       unsigned columns_number,
                                                       unsigned tuples_number,
                                                       double c_sample = 1,
                                                       double rSample = 1);

	CFDRelationData(std::unique_ptr<RelationalSchema> schema, std::vector<CFDColumnData> column_data,
                    std::vector<Transaction> data, boost::unordered_map<std::pair<int, std::string>, int, pairhash>
                    item_dict,std::vector<ItemInfo> items) :
          AbstractRelationData(std::move(schema), std::move(column_data)),
          data_rows_(std::move(data)), item_dictionary_(std::move(item_dict)),
          items_(std::move(items)) {}
};

#pragma once

#include "CSVParser.h"
#include "RelationData.h"
#include "TypedColumnData.h"

namespace model {

using TypedRelationData = AbstractRelationData<TypedColumnData>;

class ColumnLayoutTypedRelationData final : public TypedRelationData {
public:
    using TypedRelationData::AbstractRelationData;

    unsigned int GetNumRows() const final {
        if (column_data_.empty()) {
            return 0;
        } else {
            return column_data_.front().GetNumRows();
        }
    }

    static std::unique_ptr<ColumnLayoutTypedRelationData> CreateFrom(CSVParser& file_input,
                                                                     bool is_null_eq_null,
                                                                     int max_cols = -1,
                                                                     long max_rows = -1);
};

}  // namespace model

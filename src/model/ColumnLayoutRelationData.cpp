//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#include <map>
#include <memory>
#include <utility>

#include "ColumnLayoutRelationData.h"
#include "easylogging++.h"


ColumnLayoutRelationData::ColumnLayoutRelationData(std::unique_ptr<RelationalSchema> schema,
                                                   std::vector<ColumnData> columnData) :
    RelationData(std::move(schema)),
    columnData(std::move(columnData)) {}


ColumnData& ColumnLayoutRelationData::getColumnData(int columnIndex) {
    return columnData[columnIndex];
}

ColumnData const& ColumnLayoutRelationData::getColumnData(int columnIndex) const {
    return columnData[columnIndex];
}

unsigned int ColumnLayoutRelationData::getNumRows() const {
    return columnData[0].getProbingTable().size();
}

std::vector<int> ColumnLayoutRelationData::getTuple(int tupleIndex) const {
    int numColumns = schema->getNumColumns();
    std::vector<int> tuple = std::vector<int>(numColumns);
    for (int columnIndex = 0; columnIndex < numColumns; columnIndex++){
        tuple[columnIndex] = columnData[columnIndex].getProbingTableValue(tupleIndex);
    }
    return tuple;
}

/*void ColumnLayoutRelationData::shuffleColumns() {
    for (auto &columnDatum : columnData){
        columnDatum.shuffle();
    }
}*/

std::unique_ptr<ColumnLayoutRelationData> ColumnLayoutRelationData::createFrom(CSVParser &fileInput, bool isNullEqNull) {
    return createFrom(fileInput, isNullEqNull, -1, -1);
}

std::unique_ptr<ColumnLayoutRelationData> ColumnLayoutRelationData::createFrom(
        CSVParser &fileInput, bool isNullEqNull, int maxCols, long maxRows) {
    auto schema = std::make_unique<RelationalSchema>(fileInput.getRelationName(), isNullEqNull);
    std::unordered_map<std::string, int> valueDictionary;
    int nextValueId = 1;
    const int nullValueId = -1;
    int numColumns = fileInput.getNumberOfColumns();
    if (maxCols > 0) numColumns = std::min(numColumns, maxCols);
    std::vector<std::vector<int>> columnVectors = std::vector<std::vector<int>>(numColumns);
    int rowNum = 0;
    std::vector<std::string> row;

    while (fileInput.getHasNext()){
        row = fileInput.parseNext();

        if ((int)row.size() != numColumns) {
            LOG(WARNING) << "Skipping incomplete rows";
            continue;
        }

        if (maxRows <= 0 || rowNum < maxRows){
            int index = 0;
            for (std::string& field : row){
                if (field.empty()){
                    columnVectors[index].push_back(nullValueId);
                } else {
                    auto location = valueDictionary.find(field);
                    int valueId;
                    if (location == valueDictionary.end()){
                        valueDictionary[field] = nextValueId;
                        valueId = nextValueId;
                        nextValueId++;
                    } else {
                        valueId = location->second;
                    }
                    columnVectors[index].push_back(valueId);
                }
                index++;
                if (index >= numColumns) break;
            }
        } else {
            //TODO: Подумать что тут сделать
            assert(0);
        }
        rowNum++;
    }

    std::vector<ColumnData> columnData;
    for (int i = 0; i < numColumns; ++i) {
        auto column = Column(schema.get(), fileInput.getColumnName(i), i);
        schema->appendColumn(std::move(column));
        auto pli = PositionListIndex::createFor(columnVectors[i], schema->isNullEqualNull());
        columnData.emplace_back(schema->getColumn(i), std::move(pli));
    }

    schema->init();

    return std::make_unique<ColumnLayoutRelationData>(std::move(schema), std::move(columnData));
}
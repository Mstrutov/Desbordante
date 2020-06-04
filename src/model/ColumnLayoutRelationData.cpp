//
// Created by kek on 26.07.2019.
//

#include "ColumnLayoutRelationData.h"

#include <map>
#include <memory>
#include <utility>

using namespace std;

//vector <- pass by reference
ColumnLayoutRelationData::ColumnLayoutRelationData(shared_ptr<RelationalSchema>& schema, vector<shared_ptr<ColumnData>> columnData) :
    RelationData(schema),
    columnData(std::move(columnData)){}

vector<shared_ptr<ColumnData>> ColumnLayoutRelationData::getColumnData() {
    return columnData;
}

shared_ptr<ColumnData> ColumnLayoutRelationData::getColumnData(int columnIndex) {
    return columnData[columnIndex];
}

unsigned int ColumnLayoutRelationData::getNumRows() {
    return columnData[0]->getProbingTable().size();
}

vector<int> ColumnLayoutRelationData::getTuple(int tupleIndex) {
    int numColumns = schema->getNumColumns();
    vector<int> tuple = vector<int>(numColumns);
    for (int columnIndex = 0; columnIndex < numColumns; columnIndex++){
        tuple[columnIndex] = columnData[columnIndex]->getProbingTableValue(tupleIndex);
    }
    return tuple;
}

void ColumnLayoutRelationData::shuffleColumns() {
    for (auto &columnDatum : columnData){
        columnDatum->shuffle();
    }
}

shared_ptr<ColumnLayoutRelationData> ColumnLayoutRelationData::createFrom(CSVParser &fileInput, bool isNullEqNull) {
    return createFrom(fileInput, isNullEqNull, -1, -1);
}

shared_ptr<ColumnLayoutRelationData> ColumnLayoutRelationData::createFrom(CSVParser &fileInput, bool isNullEqNull, int maxCols,
                                                              long maxRows) {
    auto schema =  RelationalSchema::create(fileInput.getRelationName(), isNullEqNull);
    map<string, int> valueDictionary;
    int nextValueId = 1;
    const int nullValueId = -1;
    const int unknownValueId = 0;
    int numColumns = fileInput.getNumberOfColumns();
    if (maxCols > 0) numColumns = min(numColumns, maxCols);
    vector<vector<int>> columnVectors = vector<vector<int>>(numColumns);
    int rowNum = 0;
    vector<string> row;

    while (fileInput.getHasNext()){
        row = std::move(fileInput.parseNext());

        //skip incomplete rows
        if ((int)row.size() != numColumns) {
            continue;
        }

        if (maxRows <= 0 || rowNum < maxRows){
            int index = 0;
            for (string& field : row){
                if (field.empty()){
                    columnVectors[index].push_back(nullValueId);
                } else {
                    int valueId = valueDictionary[field];
                    if (valueId == 0){
                        valueDictionary[field] = nextValueId;
                        valueId = nextValueId;
                        nextValueId++;
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

    vector<shared_ptr<ColumnData>> columnData;
    for (int i = 0; i < numColumns; ++i) {
        auto column = make_shared<Column>(schema, fileInput.getColumnName(i), i);               //numColumns instead of i - same problem as described in Column.cpp
        schema->appendColumn(column);
        auto pli = PositionListIndex::createFor(columnVectors[i], schema->isNullEqualNull());
        auto colData = make_shared<ColumnData>(column, pli->getProbingTable(true), pli);
        columnData.emplace_back(colData);
    }
    // TODO: тут костыль: в RelationalSchema::create происходит инициализация битсета нулевого размера. Можно обойти через resize
    schema->init();

    auto ans = shared_ptr<ColumnLayoutRelationData>(new ColumnLayoutRelationData(schema, columnData));
    return ans;
}
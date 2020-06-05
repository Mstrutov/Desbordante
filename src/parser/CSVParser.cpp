#include "CSVParser.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>


using namespace std;

CSVParser::CSVParser(fs::path path): CSVParser(path, ',') {}

CSVParser::CSVParser(fs::path path, char separator) :
    source(path),
    separator(separator),
    hasNext(true),
    nextLine(),
    numberOfColumns(),
    columnNames(),
    relationName(path.filename().string()){

    //Wrong path
    if (!source) {
        std::runtime_error("Error: couldn't find file " + path.string());
    }
    // TODO: Настроить Exception
    if (separator == '\0'){
        assert(0);
    }
    getNext();
    vector<string> nextParsed = std::move(parseNext());
    numberOfColumns = nextParsed.size();
    columnNames = std::move(nextParsed);
}

/*
bool CSVParser::isSameChar(char separator, char escape) {
    return separator != '\0' && separator == escape;
}
*/

void CSVParser::getNext(){
    nextLine = "";
    getline(source, nextLine);
}

vector<string> CSVParser::parseNext() {
    vector<string> result = vector<string>();

    auto nextTokenBegin = nextLine.begin();
    auto nextTokenEnd = nextLine.begin();
    while (nextTokenEnd != nextLine.end()){
        if (*nextTokenEnd == separator){
            result.emplace_back(nextTokenBegin, nextTokenEnd);
            nextTokenBegin = nextTokenEnd + 1;
            nextTokenEnd = nextTokenBegin;
        } else {
            nextTokenEnd++;
        }
    }
    result.emplace_back(nextTokenBegin, nextTokenEnd);

    hasNext = !source.eof();
    if(hasNext){
        getNext();
    }

    return result;
}

bool CSVParser::getHasNext() { return hasNext;}
char CSVParser::getSeparator() { return separator;}
int CSVParser::getNumberOfColumns() { return numberOfColumns; }
string CSVParser::getColumnName(int index) { return columnNames[index]; }
string CSVParser::getRelationName() { return relationName; }
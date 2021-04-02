#include "CSVParser.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>


using namespace std;

inline std::string & CSVParser::rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {return !std::isspace(c);}).base(), s.end());
    return s;
}

CSVParser::CSVParser(fs::path path): CSVParser(path, ',', true) {}

CSVParser::CSVParser(fs::path path, char separator, bool hasHeader) :
    source(path),
    separator(separator),
    hasHeader(hasHeader),
    hasNext(true),
    nextLine(),
    numberOfColumns(),
    columnNames(),
    relationName(path.filename().string()){

    //Wrong path
    if (!source) {
        throw std::runtime_error("Error: couldn't find file " + path.string());
    }
    // TODO: Настроить Exception
    if (separator == '\0'){
        assert(0);
    }
    if (hasHeader) {
        getNext();
    } else {
        peekNext();
    }
    vector<string> nextParsed = std::move(parseNext());
    numberOfColumns = nextParsed.size();
    columnNames = std::move(nextParsed);
    if (!hasHeader) {
        for (int i = 0; i < numberOfColumns; i++){
            columnNames[i] = std::to_string(i);
        }
    }
}

/*
bool CSVParser::isSameChar(char separator, char escape) {
    return separator != '\0' && separator == escape;
}
*/

void CSVParser::getNext(){
    nextLine = "";
    getline(source, nextLine);
    rtrim(nextLine);
}

void CSVParser::peekNext(){
    int len = source.tellg();
    getNext();
    source.seekg(len, std::ios_base::beg);
}

vector<string> CSVParser::parseNext() {
    vector<string> result = vector<string>();

    auto nextTokenBegin = nextLine.begin();
    auto nextTokenEnd = nextLine.begin();
    bool isEscaped = false;
    while (nextTokenEnd != nextLine.end()){
        if (!isEscaped && *nextTokenEnd == separator){
            result.emplace_back(nextTokenBegin, nextTokenEnd);
            nextTokenBegin = nextTokenEnd + 1;
            nextTokenEnd = nextTokenBegin;
        } else {
            isEscaped ^= (*nextTokenEnd == escapeSymbol);
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
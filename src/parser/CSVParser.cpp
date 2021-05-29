#include "CSVParser.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

inline std::string & CSVParser::rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {return !std::isspace(c);}).base(), s.end());
    return s;
}

CSVParser::CSVParser(const std::filesystem::path& path): CSVParser(path, ',', true) {}

CSVParser::CSVParser(const std::filesystem::path& path, char separator, bool hasHeader) :
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
    std::vector<std::string> nextParsed = parseNext();
    numberOfColumns = nextParsed.size();
    columnNames = std::move(nextParsed);
    if (!hasHeader) {
        for (int i = 0; i < numberOfColumns; i++){
            columnNames[i] = std::to_string(i);
        }
    }
}

void CSVParser::getNext(){
    nextLine = "";
    getline(source, nextLine);
    rtrim(nextLine);
}

void CSVParser::peekNext() {
    int len = source.tellg();
    getNext();
    source.seekg(len, std::ios_base::beg);
}

std::vector<std::string> CSVParser::parseNext() {
    std::vector<std::string> result = std::vector<std::string>();

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
    if (nextTokenBegin != nextLine.begin() || nextTokenBegin != nextTokenEnd) {
        result.emplace_back(nextTokenBegin, nextTokenEnd);
    }

    hasNext = !source.eof();
    if(hasNext){
        getNext();
    }

    return result;
}
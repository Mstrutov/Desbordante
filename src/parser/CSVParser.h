//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using std::vector, std::string;

class CSVParser {
private:
    std::ifstream source;
    char separator;
    char escapeSymbol = '\"';
    bool hasHeader;
    bool hasNext;
    std::string nextLine;
    int numberOfColumns;
    vector<string> columnNames;
    string relationName;
    void getNext();
    void peekNext();

    static inline std::string & rtrim(std::string &s);
public:
    CSVParser(fs::path path);
    CSVParser(fs::path path, char separator, bool hasHeader);
    //bool isSameChar(char separator, char escape);
    vector<string> parseNext();
    bool getHasNext();
    char getSeparator();
    int getNumberOfColumns();
    string getColumnName(int index);
    string getRelationName();
};
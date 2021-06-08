//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

class CSVParser {
private:
    std::ifstream source;
    char separator;
    char escapeSymbol = '\"';
    bool hasHeader;
    bool hasNext;
    std::string nextLine;
    int numberOfColumns;
    std::vector<std::string> columnNames;
    std::string relationName;
    void getNext();
    void peekNext();

    static inline std::string & rtrim(std::string &s);
public:
    explicit CSVParser(const std::filesystem::path& path);
    CSVParser(const std::filesystem::path& path, char separator, bool hasHeader);
    //bool isSameChar(char separator, char escape);
    std::vector<std::string> parseNext();
    bool getHasNext() const { return hasNext;}
    char getSeparator() const { return separator;}
    int getNumberOfColumns() const { return numberOfColumns; }
    std::string getColumnName(int index) const { return columnNames[index]; }
    std::string getRelationName() const { return relationName; }
};

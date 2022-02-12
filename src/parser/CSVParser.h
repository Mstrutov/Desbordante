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
    std::ifstream source_;
    char separator_;
    char escape_symbol_ = '\"';
    bool has_header_;
    bool has_next_;
    std::string next_line_;
    int number_of_columns_;
    std::vector<std::string> column_names_;
    std::string relation_name_;
    void GetNext();
    void PeekNext();

    static inline std::string& rtrim(std::string& s);

public:
    CSVParser() = default;
    explicit CSVParser(const std::filesystem::path& path);
    CSVParser(const std::filesystem::path& path, char separator, bool has_header);
    //bool isSameChar(char separator, char escape);
    std::vector<std::string> ParseNext();
    bool GetHasNext() const {
        return has_next_;
    }
    char GetSeparator() const {
        return separator_;
    }
    int GetNumberOfColumns() const { return number_of_columns_; }
    std::string GetColumnName(int index) const { return column_names_[index]; }
    std::string GetRelationName() const { return relation_name_; }
};

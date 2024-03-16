//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "model/table/idataset_stream.h"

struct CSVConfig {
    std::filesystem::path path;
    char separator;
    bool has_header;
};

class CSVParser : public model::IDatasetStream {
private:
    std::ifstream source_;
    char separator_;
    char escape_symbol_ = '\\';
    char quote_ = '\"';
    bool has_header_;
    bool has_next_;
    std::string next_line_;
    int number_of_columns_;
    std::vector<std::string> column_names_;
    std::string relation_name_;
    void GetNext();
    void PeekNext();
    void GetLine(unsigned long long const line_index);
    std::vector<std::string> ParseString(std::string const& s) const;
    void GetNextIfHas();
    void SkipLine();
    std::optional<char> DeduceSeparator();
    bool CheckSeparator(char sep);

    inline static std::string& rtrim(std::string& s);

public:
    CSVParser() = default;
    explicit CSVParser(std::filesystem::path const& path);
    CSVParser(std::filesystem::path const& path, char separator, bool has_header);
    explicit CSVParser(CSVConfig const& csv_config);

    std::vector<std::string> GetNextRow() override;
    std::string GetUnparsedLine(unsigned long long const line_index);
    std::vector<std::string> ParseLine(unsigned long long const line_index);

    std::optional<char> ValidateSeparator();

    bool HasNextRow() const override {
        return has_next_;
    }

    char GetSeparator() const {
        return separator_;
    }

    size_t GetNumberOfColumns() const override {
        return number_of_columns_;
    }

    std::string GetColumnName(size_t index) const override {
        return column_names_[index];
    }

    std::string GetRelationName() const override {
        return relation_name_;
    }

    void Reset() override;
};

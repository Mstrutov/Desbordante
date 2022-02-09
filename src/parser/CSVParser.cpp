#include "CSVParser.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

inline std::string& CSVParser::rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) { return !std::isspace(c); }).base(),
            s.end());
    return s;
}

CSVParser::CSVParser(const std::filesystem::path& path) : CSVParser(path, ',', true) {}

CSVParser::CSVParser(const std::filesystem::path& path, char separator, bool has_header)
    : source_(path),
      separator_(separator),
      has_header_(has_header),
      has_next_(true),
      next_line_(),
      number_of_columns_(),
      column_names_(),
      relation_name_(path.filename().string()) {
    //Wrong path
    if (!source_) {
        throw std::runtime_error("Error: couldn't find file " + path.string());
    }
    // TODO: Настроить Exception
    if (separator == '\0') {
        assert(0);
    }
    if (has_header) {
        GetNext();
    } else {
        PeekNext();
    }
    std::vector<std::string> next_parsed = ParseNext();
    number_of_columns_ = next_parsed.size();
    column_names_ = std::move(next_parsed);
    if (!has_header) {
        for (int i = 0; i < number_of_columns_; i++) {
            column_names_[i] = std::to_string(i);
        }
    }
}

void CSVParser::GetNext() {
    next_line_ = "";
    getline(source_, next_line_);
    rtrim(next_line_);
}

void CSVParser::PeekNext() {
    int len = source_.tellg();
    GetNext();
    source_.seekg(len, std::ios_base::beg);
}

std::vector<std::string> CSVParser::ParseNext() {
    std::vector<std::string> result = std::vector<std::string>();

    auto next_token_begin = next_line_.begin();
    auto next_token_end = next_line_.begin();
    bool is_escaped = false;
    while (next_token_end != next_line_.end()) {
        if (!is_escaped && *next_token_end == separator_) {
            result.emplace_back(next_token_begin, next_token_end);
            next_token_begin = next_token_end + 1;
            next_token_end = next_token_begin;
        } else {
            is_escaped ^= (*next_token_end == escape_symbol_);
            next_token_end++;
        }
    }
    if (next_token_begin != next_line_.begin() || next_token_begin != next_token_end) {
        result.emplace_back(next_token_begin, next_token_end);
    }

    has_next_ = !source_.eof();
    if (has_next_) {
        GetNext();
    }

    return result;
}
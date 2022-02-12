#include "CSVParser.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

inline std::string& CSVParser::rtrim(std::string& s) {
    boost::trim_right(s);
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
    std::getline(source_, next_line_);
    rtrim(next_line_);
}

void CSVParser::PeekNext() {
    int len = source_.tellg();
    GetNext();
    source_.seekg(len, std::ios_base::beg);
}

void CSVParser::GetLine(const unsigned long long line_index) {
    source_.clear();
    source_.seekg(0);

    next_line_ = "";

    /* Skip header */
    if (has_header_) {
        std::getline(source_, next_line_);
    }

    /* Index is less than the line number by one */
    for (unsigned long long i = 0; i < line_index + 1; ++i) {
        std::getline(source_, next_line_);
    }

    rtrim(next_line_);
}

void CSVParser::GetNextIfHas() {
    has_next_ = !source_.eof();
    if (has_next_) {
        GetNext();
    }
}

std::string CSVParser::GetUnparsedLine(const unsigned long long line_index) {
    GetLine(line_index);
    std::string line = next_line_;

    /* For correct work of ParseNext() after this method */
    GetNextIfHas();

    return line;
}

std::vector<std::string> CSVParser::ParseString(const std::string& s) {
    std::vector<std::string> tokens;
    boost::escaped_list_separator<char> list_sep(escape_symbol_, separator_, quote_);
    boost::tokenizer<boost::escaped_list_separator<char>> tokenizer(s, list_sep);

    for (auto& token : tokenizer) {
        tokens.push_back(std::move(token));
    }

    return tokens;
}

std::vector<std::string> CSVParser::ParseLine(const unsigned long long line_index) {
    std::string unparsed_line = GetUnparsedLine(line_index);
    std::vector<std::string> parsed = ParseString(unparsed_line);
    return parsed;
}

std::vector<std::string> CSVParser::ParseNext() {
    std::vector<std::string> result = ParseString(next_line_);

    GetNextIfHas();

    return result;
}

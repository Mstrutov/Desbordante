#include "csv_parser.h"

#include <cassert>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

inline std::string& CSVParser::Rtrim(std::string& s) {
    if (!s.empty() && s.back() == '\n') s.pop_back();
    return s;
}

CSVParser::CSVParser(std::filesystem::path const& path) : CSVParser(path, ',', true) {}

CSVParser::CSVParser(std::filesystem::path const& path, char separator, bool has_header)
    : source_(path),
      separator_(separator),
      has_header_(has_header),
      has_next_(true),
      next_line_(),
      number_of_columns_(),
      column_names_(),
      relation_name_(path.filename().string()) {
    // Wrong path
    if (!source_) {
        throw std::runtime_error("Error: couldn't find file " + path.string());
    }
    if (separator == '\0') {
        throw std::invalid_argument("Invalid separator");
    }
    if (has_header) {
        GetNext();
    } else {
        PeekNext();
    }

    std::vector<std::string> next_parsed = CSVParser::GetNextRow();
    number_of_columns_ = next_parsed.size();
    column_names_ = std::move(next_parsed);

    if (!has_header) {
        for (int i = 0; i < number_of_columns_; i++) {
            column_names_[i] = std::to_string(i);
        }
    }
}

CSVParser::CSVParser(CSVConfig const& csv_config)
    : CSVParser(csv_config.path, csv_config.separator, csv_config.has_header) {}

void CSVParser::GetNext() {
    next_line_ = "";
    std::getline(source_, next_line_);
    Rtrim(next_line_);
}

void CSVParser::PeekNext() {
    int len = source_.tellg();
    GetNext();
    source_.seekg(len, std::ios_base::beg);
}

void CSVParser::SkipLine() {
    source_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void CSVParser::Reset() {
    source_.clear();
    source_.seekg(0);

    next_line_.clear();
    has_next_ = true;

    // Skip header
    if (has_header_) {
        SkipLine();
    }

    // For correctness of GetNextRow() after this method
    GetNextIfHas();
}

void CSVParser::GetLine(unsigned long long const line_index) {
    Reset();

    // Index is less than the line number by one. Skip line_index lines
    for (unsigned long long i = 0; i < line_index; ++i) {
        SkipLine();
    }

    std::getline(source_, next_line_);

    Rtrim(next_line_);
}

void CSVParser::GetNextIfHas() {
    has_next_ = !source_.eof();

    if (has_next_) {
        if (source_.peek() == std::ifstream::traits_type::eof()) {  // Check for the last newline
            has_next_ = false;
            return;
        }
        GetNext();
    }
}

std::string CSVParser::GetUnparsedLine(unsigned long long const line_index) {
    GetLine(line_index);
    std::string line = next_line_;

    // For correctness of GetNextRow() after this method
    GetNextIfHas();

    return line;
}

std::vector<std::string> CSVParser::ParseString(std::string const& s) const {
    std::size_t const length = s.size();
    std::string t;
    for (std::size_t index = 0; index < length; ++index) {
        if (s[index] == '\\') {  // preserve \ through boost parsing
            t.append(2, '\\');
        } else if (s[index] == '"') {  // preserve " through boost parsing
            t.append("\\\"\"");
        } else {
            t.push_back(s[index]);
        }
    }

    std::vector<std::string> tokens;
    tokens.reserve(number_of_columns_);
    boost::escaped_list_separator<char> list_sep(escape_symbol_, separator_, quote_);
    boost::tokenizer<boost::escaped_list_separator<char>> tokenizer(t, list_sep);

    for (auto& token : tokenizer) {
        std::size_t const token_length = token.size();
        bool is_enclosed =
                token_length >= 2 && token.front() == '"' &&
                token.back() == '"';  // states whether a field is enclosed in double-quotes

        std::string new_token;
        for (std::size_t index = 0; index < token_length; ++index) {
            if (token[index] == '"') {
                if (is_enclosed && index > 0 && index < token_length - 2 &&
                    token[index + 1] == '"') {  // transfer "" to " if the current field is enclosed
                                                // in double-quotes
                    new_token.push_back(token[index]);
                    ++index;
                }
            } else {
                new_token.push_back(token[index]);
            }
        }
        tokens.push_back(std::move(new_token));
    }

    return tokens;
}

std::vector<std::string> CSVParser::ParseLine(unsigned long long const line_index) {
    std::string unparsed_line = GetUnparsedLine(line_index);
    std::vector<std::string> parsed = ParseString(unparsed_line);
    return parsed;
}

std::vector<std::string> CSVParser::GetNextRow() {
    std::vector<std::string> result = ParseString(next_line_);
    if (number_of_columns_ == 1 && result.empty()) {
        result = {""};
    }

    GetNextIfHas();

    return result;
}

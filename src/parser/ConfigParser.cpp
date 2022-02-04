//
// Created by kek on 10.07.19.
//

#include "ConfigParser.h"

#include <fstream>

#include "json.hpp"

using namespace std;
using json = nlohmann::json;

ConfigParser::ConfigParser(fs::path path) {
    ifstream stream(path);
    json json = json::parse(string(istreambuf_iterator<char>(stream), istreambuf_iterator<char>()));
    input_path_ = json["input_path"];
    output_path_ = json["output_path"];
}

string ConfigParser::GetOutputPath() {
    return output_path_;
}

string ConfigParser::GetInputPath() {
    return input_path_;
}
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
    inputPath = json["inputPath"];
    outputPath = json["outputPath"];
}

string ConfigParser::getOutputPath() {
    return outputPath;
}

string ConfigParser::getInputPath() {
    return inputPath;
}
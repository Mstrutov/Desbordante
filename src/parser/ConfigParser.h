//
// Created by kek on 10.07.19.
//

#pragma once

//now its not experimental - just "filesystem". Maybe ifndef - experimental or not experimental, if .
//
#include <experimental/filesystem>
#include <string>

namespace fs = std::experimental::filesystem;


class ConfigParser {
private:
    std::string outputPath;
    std::string inputPath;

public:
    ConfigParser(fs::path path);
    std::string getOutputPath();
    std::string getInputPath();

};
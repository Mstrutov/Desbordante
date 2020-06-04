//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

// macro: filesystem - defined ? filesystem : experimental/filesystem
#include <filesystem>
#include <string>

namespace fs = std::filesystem;


class ConfigParser {
private:
    std::string outputPath;
    std::string inputPath;

public:
    ConfigParser(fs::path path);
    std::string getOutputPath();
    std::string getInputPath();

};
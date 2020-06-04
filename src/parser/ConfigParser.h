//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

//now its not experimental - just "filesystem". Maybe ifndef - experimental or not experimental, if .
//
#ifdef #include <filesystem>
#define
#endif
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
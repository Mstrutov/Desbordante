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
    std::string output_path_;
    std::string input_path_;

public:
    ConfigParser(fs::path path);
    std::string GetOutputPath();
    std::string GetInputPath();

};
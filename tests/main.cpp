//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#include <filesystem>
#include <string>

#include "logging/easylogging++.h"
#include "gtest/gtest.h"

INITIALIZE_EASYLOGGINGPP

std::string GetSelfpath();

int main(int argc, char **argv) {
    el::Configurations default_conf;
    default_conf.setToDefault();
    default_conf.set(el::Level::Global,
                     el::ConfigurationType::Enabled, "false");
    el::Loggers::reconfigureAllLoggers(default_conf);

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

//Overall, just give file path to csv to the main when execute via console or smth.
//https://stackoverflow.com/a/27460370 - better use this function
std::string GetSelfpath(){
    /*std::vector<char> buf(400);
    unsigned int len;

    do {
        buf.resize(buf.size() + 100);
        // подстава: readlink не работает для шиндовс
        len = ::readlink("/proc/self/exe", &(buf[0]), buf.size());
    } while (buf.size() == len);

    if (len > 0)
    {
        buf[len] = '\0';
        return (std::string(&(buf[0])));

    */
    return std::filesystem::current_path().string();
}

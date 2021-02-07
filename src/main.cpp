#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>

#include <boost/program_options.hpp>

#include "logging/easylogging++.h"

#include "ColumnLayoutRelationData.h"
#include "ConfigParser.h"
#include "algorithms/Pyro.h"
#include "algorithms/TaneX.h"

namespace po = boost::program_options;

INITIALIZE_EASYLOGGINGPP

int main(int argc, char const *argv[]) {
    std::string alg;
    std::string dataset;
    char separator = ',';
    bool hasHeader = true;
    int seed = 0;
    double error = 0.01;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "print help")
        ("algo", po::value<string>(&alg), "algorithm [pyro|tane]")
        ("data", po::value<string>(&dataset), "path to dataset CSV file")
        ("sep", po::value<char>(&separator), "CSV separator")
        ("hasHeader", po::value<bool>(&hasHeader), "CSV header presence flag [true|false]. Default true")
        ("seed", po::value<int>(&seed), "RNG seed")
        ("error", po::value<double>(&error), "error for AFD algorithms. Default 0.01")
    ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (po::error &e) {
        std::cout << e.what() << std::endl;
        return 0;
    }

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    el::Loggers::configureFromGlobal("logging.conf");

    std::transform(alg.begin(), alg.end(), alg.begin(), [](unsigned char c){ return std::tolower(c); });
    std::cout << "Input: algorithm \"" << alg
              << "\" with seed " << std::to_string(seed)
              << ", error \"" << std::to_string(error)
              << "\" and dataset \"" << dataset
              << "\" with separator \'" << separator
              << "\'. Header is " << (hasHeader ? "" : "not ") << "present. " << std::endl;
    auto path = std::filesystem::current_path() / "inputData" / dataset;
    if (alg == "pyro") {
        try {
            Pyro algInstance(path, separator, hasHeader, seed, error);
            double elapsedTime = algInstance.execute();
            std::cout << "> ELAPSED TIME: " << elapsedTime << std::endl;
        } catch (std::runtime_error& e) {
            std::cout << e.what() << std::endl;
            return 1;
        }
    } else if (alg == "tane"){
        try {
            Tane algInstance(path, separator, hasHeader);
            algInstance.execute();
        } catch (std::runtime_error& e) {
            std::cout << e.what() << std::endl;
            return 1;
        }
    } else {
        std::cout << "Error - no matching algorithm. Available algorithms are:\n\tpyro\n\ttane\n" << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }
    return 0;
}

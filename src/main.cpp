#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <easylogging++.h>

#include "algo_factory.h"

namespace po = boost::program_options;

INITIALIZE_EASYLOGGINGPP

template<typename BetterEnumType>
static std::string EnumToAvailableValues() {
    std::stringstream avail_values;

    avail_values << '[';

    for (auto const& name : BetterEnumType::_names()) {
        avail_values << name << '|';
    }

    avail_values.seekp(-1, avail_values.cur);
    avail_values << ']';

    return avail_values.str();
}

static bool CheckOptions(std::string const& alg, double error) {
    if (!algos::Algo::_is_valid(alg.c_str())) {
        std::cout << "ERROR: no matching algorithm."
                     " Available algorithms are:\n" +
                     EnumToAvailableValues<algos::Algo>() + '\n';
        return false;
    }

    if (error > 1 || error < 0) {
        std::cout << "ERROR: error should be between 0 and 1.\n";
    }

    return true;
}

int main(int argc, char const* argv[]) {
    std::string alg;
    std::string dataset;
    char separator = ',';
    bool has_header = true;
    int seed = 0;
    double error = 0.0;
    unsigned int max_lhs = -1;
    unsigned int parallelism = 0;

    std::string const algo_desc = "algorithm to use. Available algorithm:\n" +
                                  EnumToAvailableValues<algos::Algo>() +
                                  " for FD mining.";

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "print this message and exit")
        ("mine", "type of dependency to mine")
        ("algo", po::value<std::string>(&alg),
         algo_desc.c_str())
        ("data", po::value<std::string>(&dataset),
         "path to CSV file, relative to ./inputData")
        ("sep", po::value<char>(&separator)->default_value(separator),
         "CSV separator")
        ("hasHeader", po::value<bool>(&has_header)->default_value(has_header),
         "CSV header presence flag [true|false]")
        ("seed", po::value<int>(&seed), "RNG seed")
        ("error", po::value<double>(&error)->default_value(error),
         "error for AFD algorithms")
        ("maxLHS", po::value<unsigned int>(&max_lhs)->default_value(max_lhs),
         "max considered LHS size")
        ("threads", po::value<unsigned int>(&parallelism),
         "number of threads to use")
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

    std::transform(alg.begin(), alg.end(), alg.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (!CheckOptions(alg, error)) {
        std::cout << desc << std::endl;
        return 1;
    }

    std::cout << "Input: algorithm \"" << alg
              << "\" with seed " << std::to_string(seed)
              << ", error \"" << std::to_string(error)
              << ", maxLHS \"" << std::to_string(max_lhs)
              << "\" and dataset \"" << dataset
              << "\" with separator \'" << separator
              << "\'. Header is " << (has_header ? "" : "not ") << "present. " << std::endl;



    std::unique_ptr<FDAlgorithm> algorithm_instance =
        algos::CreateAlgorithmInstance(alg, dataset, separator, has_header, error, max_lhs,
                                       parallelism, seed);

    try {
        unsigned long long elapsed_time = algorithm_instance->Execute();
        std::cout << "> ELAPSED TIME: " << elapsed_time << std::endl;
    } catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}

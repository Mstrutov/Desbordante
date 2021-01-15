//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>

#include "logging/easylogging++.h"

#include "ColumnLayoutRelationData.h"
#include "ConfigParser.h"
#include "algorithms/Pyro.h"
#include "algorithms/TaneX.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char const *argv[]) {
    el::Loggers::configureFromGlobal("logging.conf");

    using std::cout, std::endl;
    std::vector<std::string> strs(argv + 1, argv + argc);

    std::string alg;
    std::string dataset;
    char separator = ',';
    bool hasHeader = true;
    int seed = 0;
    std::vector<std::string> availableAlgs {"pyro", "tane"};
    for (auto& str : strs) {
        if (str.size() > 2) {
            if (str.find("--algo=") == 0) {
                alg = str.substr(7);
            } else if (str.find("--data=") == 0) {
                dataset = str.substr(7);
            } else if (str.find("--sep=") == 0) {
                separator = str[6];
            } else if (str.find("--hasHeader=") == 0) {
                hasHeader = (std::tolower(str[12]) == 't');
            } else if (str.find("--seed=") == 0) {
                seed = std::stoi(str.substr(7));
            }
        }
    }
    transform(alg.begin(), alg.end(), alg.begin(), [](unsigned char c){ return std::tolower(c); });
    cout << "Input: algorithm \"" << alg
         << "\" with seed " << std::to_string(seed)
         << " and dataset \"" << dataset
         << "\" with separator \'" << separator
         << "\'. Header is " << (hasHeader ? "" : "not ") << "present. " << endl;
    auto path = std::filesystem::current_path() / "inputData" / dataset;
    if (alg.empty() || dataset.empty()) {
        cout << "Couldn't recognize the algorithm and the dataset.\n" <<
                "Try launching with options: --algo=<name>    --data=<../relative/path/to/dataset>\n" <<
                "                            --sep='<symbol>' --hasHeader=[t|f|T|F] \n" <<
                "                            --seed='<int>' \n" <<
                "Available algorithms are:\n\tpyro\n\ttane\n";
    } else if (alg == "pyro") {
        try {
            Pyro algInstance(path, separator, hasHeader, seed);
            double elapsedTime = algInstance.execute();
            cout << "> ELAPSED TIME: " << elapsedTime << endl;
        } catch (std::runtime_error& e) {
            cout << e.what() << endl;
            return 1;
        }
    } else if (alg == "tane"){
        try {
            Tane algInstance(path, separator, hasHeader);
            algInstance.execute();
        } catch (std::runtime_error& e) {
            cout << e.what() << endl;
            return 1;
        }
    } else {
        cout << "Error - no matching algorithm. Available algorithms are:\n\tpyro\n\ttane\n" << endl;
        return 1;
    }
    return 0;
}

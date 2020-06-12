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

#include "ColumnLayoutRelationData.h"
#include "ConfigParser.h"
#include "algorithms/Pyro.h"
#include "algorithms/TaneX.h"

int main(int argc, char const *argv[]) {
    using std::cout, std::endl;
    std::vector<std::string> strs(argv + 1, argv + argc);

    std::string alg;
    std::string dataset;
    std::vector<std::string> availableAlgs {"pyro", "tane"};
    for (auto& str : strs) {
        if (str.size() > 2) {
            if (str.find("-algo=") == 0) {
                alg = str.substr(6);
            } else if (str.find("-data=") == 0) {
                dataset = str.substr(6);
            }
        }
    }
    transform(alg.begin(), alg.end(), alg.begin(), [](unsigned char c){ return std::tolower(c); });
    cout << "Input: algorithm \"" << alg  << "\" and dataset \"" << dataset << "\"" << endl;
    auto path = std::filesystem::current_path() / "inputData" / dataset;
    if (alg.empty() || dataset.empty()) {
        cout << "Couldn't recognize the algorithm and the dataset.\n" <<
                "Try launching with options: -algo=<name> -data=<../relative/path/to/dataset>\n" <<
                "Available algorithms are:\n\tpyro\n\ttane\n";
    } else if (alg == "pyro") {
        Pyro algInstance(path);
        try {
            algInstance.execute();
        } catch (std::runtime_error& e) {
            cout << e.what() << endl;
            return 1;
        }
    } else if (alg == "tane"){
        Tane algInstance(path);
        try {
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

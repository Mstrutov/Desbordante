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
    for (auto& str : strs) {
        cout << str << endl;
    }
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
    auto path = std::filesystem::current_path() / dataset;
    if (alg == "pyro") {
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
        cout << "Error: available algorithms are: - pyro\n -tane\n" << endl;
        return 1;
    }
    return 0;
}

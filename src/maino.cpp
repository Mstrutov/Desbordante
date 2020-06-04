//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#include "ColumnLayoutRelationData.h"
#include "ConfigParser.h"

using namespace std;

int maino(){
    ConfigParser configParser("config.json");
    CSVParser parser(configParser.getInputPath());
    return 0;
}

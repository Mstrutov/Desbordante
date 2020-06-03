//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#include "model/ColumnLayoutRelationData.h"
#include "parser/ConfigParser.h"

using namespace std;

int maino(){
    ConfigParser configParser("config.json");
    CSVParser parser(configParser.getInputPath());
    return 0;
}

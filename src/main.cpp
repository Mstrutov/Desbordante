#include "model/ColumnLayoutRelationData.h"
#include "parser/ConfigParser.h"

using namespace std;

int main(){
    ConfigParser configParser("config.json");
    CSVParser parser(configParser.getInputPath());
    return 0;
}

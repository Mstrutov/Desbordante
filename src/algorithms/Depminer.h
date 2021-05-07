#pragma once

#include "CSVParser.h"
#include "FDAlgorithm.h"

class Depminer : public FDAlgorithm {
public:
    
    explicit Depminer(std::filesystem::path const& path, char separator = ',', bool hasHeader = true) 
        : FDAlgorithm(path, separator, hasHeader) {}

    unsigned long long execute() override;

};

/*

Что написано в статье 
"Efficient Discovery..." Jean-Mart petit, Lotfi Lakhal

Striped partition database
AgreeSets
Maximal sets
Complement sets
LHS of FDs from Maximal sets



*/
//когда начинать считать время в execute? до получения базы или после

#pragma once

#include <vector>
#include <string>

struct Dataset {
    std::string name;
    unsigned int hash;
    char separator;
    bool header_presence;
};

class DatasetCollection {
protected:
    virtual size_t DatasetQuantity() = 0;
    virtual std::string DatasetName(size_t i) = 0;
    virtual char Separator(size_t i) = 0;
    virtual bool HasHeader(size_t i) = 0;
    virtual unsigned int Hash(size_t i) = 0;
};

class LightDatasets : public DatasetCollection {
protected:
    std::vector<Dataset> datasets_ = {
            {"CIPublicHighway10k.csv", 21537, ',', true},
            {"neighbors10k.csv", 62210, ',', true},
            {"WDC_astronomical.csv", 49417, ',', true},
            {"WDC_age.csv", 44913, ',', true},
            {"WDC_appearances.csv", 34497, ',', true},
            {"WDC_astrology.csv", 12109, ',', true},
            {"WDC_game.csv", 32103, ',', true},
            {"WDC_science.csv", 44913, ',', true},
            {"WDC_symbols.csv", 2160, ',', true},
            {"breast_cancer.csv", 58641, ',', true},
            {"WDC_kepler.csv", 6284, ',', true}
    };

    // DEPRECATED -- just use
    // for (auto dataset : LightDatasets::datasets) { ... }
    size_t DatasetQuantity() override { return datasets_.size(); }
    std::string DatasetName(size_t i) override { return datasets_[i].name; }
    char Separator(size_t i) override { return datasets_[i].separator; }
    bool HasHeader(size_t i) override { return datasets_[i].header_presence; }
    unsigned int Hash(size_t i) override { return datasets_[i].hash; }
};

class HeavyDatasets : public DatasetCollection {
protected:
    std::vector<Dataset> datasets_ = {
            {"adult.csv", 20873, ';', false},
            {"CIPublicHighway.csv", 32696, ',', true},
            {"EpicMeds.csv", 19617, '|', true},
            {"EpicVitals.csv", 24082, '|', true},
            {"iowa1kk.csv", 55285, ',', true},
            {"LegacyPayors.csv", 34874, '|', true}
    };

    // DEPRECATED -- just use
    // for (auto dataset : HeavyDatasets::datasets) { ... }
    size_t DatasetQuantity() override { return datasets_.size(); }
    std::string DatasetName(size_t i) override { return datasets_[i].name; }
    char Separator(size_t i) override { return datasets_[i].separator; }
    bool HasHeader(size_t i) override { return datasets_[i].header_presence; }
    unsigned int Hash(size_t i) override { return datasets_[i].hash; }
};
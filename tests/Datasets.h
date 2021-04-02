#include <vector>
#include <string>

class DatasetCollection {
protected:
    virtual size_t datasetQuantity() = 0;
    virtual std::string dataset(size_t i) = 0;
    virtual char separator(size_t i) = 0;
    virtual bool hasHeader(size_t i) = 0;
    virtual unsigned int hash(size_t i) = 0;
};

class LightDatasets : public DatasetCollection {
protected:
    std::vector<unsigned int> hashes = {43667, 62210, 49417, 44913, 34497, 12109, 32103, 44913, 2160, 6284};
    std::vector<std::string> datasets = {"CIPublicHighway50k.csv",
                               "neighbors100k.csv",
                               "WDC_astronomical.csv",
                               "WDC_age.csv",
                               "WDC_appearances.csv",
                               "WDC_astrology.csv",
                               "WDC_game.csv",
                               "WDC_science.csv",
                               "WDC_symbols.csv",
                               "WDC_kepler.csv"};
    std::vector<char> separators = {',', ',', ',', ',', ',', ',', ',', ',', ',', ','};
    std::vector<bool> header_presence = {true, true, true, true, true, true, true, true, true, true};

    size_t datasetQuantity() override { return datasets.size(); }
    std::string dataset(size_t i) override { return datasets[i]; }
    char separator(size_t i) override { return separators[i]; }
    bool hasHeader(size_t i) override { return header_presence[i]; }
    unsigned int hash(size_t i) override { return hashes[i]; }
};

class HeavyDatasets : public DatasetCollection {
protected:
    std::vector<unsigned int> hashes = {20873, 58641, 32696, 19617, 24082, 55285, 34874, 62210};
    std::vector<std::string> datasets = {"adult.csv", "breast_cancer.csv", "CIPublicHighway.csv", "EpicMeds.csv",
                               "EpicVitals.csv", "iowa1kk.csv", "LegacyPayors.csv", "neighbors100k.csv"};
    std::vector<char> separators = {';', ',', ',', '|', '|', ',', '|', ','};
    std::vector<bool> headerPresence = {false, true, true, true, true, true, true, true};

    size_t datasetQuantity() override { return datasets.size(); }
    std::string dataset(size_t i) override { return datasets[i]; }
    char separator(size_t i) override { return separators[i]; }
    bool hasHeader(size_t i) override { return headerPresence[i]; }
    unsigned int hash(size_t i) override { return hashes[i]; }
};
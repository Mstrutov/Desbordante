#include <chrono>

#include "FdG1Strategy.h"
#include "KeyG1Strategy.h"
#include "Pyro.h"

void Pyro::execute() {
    auto relation = ColumnLayoutRelationData::createFrom(inputGenerator_, configuration_.isNullEqualNull);
    auto schema = relation->getSchema();
    auto profilingContext = std::make_shared<ProfilingContext>(
            configuration_,
            relation,
            uccConsumer_,
            fdConsumer_,
            cachingMethod_,
            evictionMethod_,
            cachingMethodValue
            );

    auto startTime = std::chrono::system_clock::now();

    std::function<bool(DependencyCandidate const&, DependencyCandidate const&)> launchPadOrder;
    if (configuration_.launchPadOrder == "arity") {
        launchPadOrder = DependencyCandidate::fullArityErrorComparator;
    } else if (configuration_.launchPadOrder == "error") {
        launchPadOrder = DependencyCandidate::fullErrorArityComparator;
    } else {
        throw std::runtime_error("Unknown comparator type");
    }

    int nextId = 0;
    if (configuration_.isFindKeys) {
        std::shared_ptr<DependencyStrategy> strategy;
        if (configuration_.uccErrorMeasure == "g1prime") {
            strategy = std::dynamic_pointer_cast<DependencyStrategy>(
                    std::make_shared<KeyG1Strategy>(configuration_.maxUccError, configuration_.errorDev));
        } else {
            throw std::runtime_error("Unknown key error measure.");
        }
        searchSpaces_.push_back(std::make_shared<SearchSpace>(nextId++, strategy, schema, launchPadOrder));
    }
    if (configuration_.isFindFds) {
        for (auto rhs : schema->getColumns()) {
            std::shared_ptr<DependencyStrategy> strategy;
            if (configuration_.uccErrorMeasure == "g1prime") {
                strategy = std::dynamic_pointer_cast<DependencyStrategy>(
                        std::make_shared<FdG1Strategy>(rhs, configuration_.maxUccError, configuration_.errorDev));
            } else {
                throw std::runtime_error("Unknown key error measure.");
            }
            searchSpaces_.push_back(std::make_shared<SearchSpace>(nextId++, strategy, schema, launchPadOrder));
        }
    }
    unsigned int totalErrorCalcCount = 0;
    unsigned long long totalAscension = 0;
    unsigned long long totalTrickle = 0;
    for (auto& searchSpace : searchSpaces_) {
        searchSpace->setContext(profilingContext);
        searchSpace->ensureInitialized();
        searchSpace->discover();
        searchSpace->printStats();
        totalErrorCalcCount += searchSpace->getErrorCalcCount();
        totalAscension += searchSpace->ascending / 1000000;
        totalTrickle += searchSpace->tricklingDown/ 1000000;
    }
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);

    std::cout << "Time: " << elapsed_milliseconds.count() << " milliseconds" << std::endl;
    std::cout << "Error calculation count: " << totalErrorCalcCount << std::endl;
    std::cout << "Total ascension: " << totalAscension << std::endl;
    std::cout << "Total trickle: " << totalTrickle << std::endl;

}

Pyro::Pyro(fs::path const &path) : inputGenerator_(path) {
    uccConsumer_ = [](auto const& key) { std::cout << "Found ucc: " << key.toString() << std::endl; };
    fdConsumer_ = [](auto const& key) { std::cout << "Found fd: " << key.toString() << std::endl; };
}

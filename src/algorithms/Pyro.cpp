#include <chrono>

#include "FdG1Strategy.h"
#include "KeyG1Strategy.h"
#include "Pyro.h"

double Pyro::execute() {

    auto relation = ColumnLayoutRelationData::createFrom(inputGenerator_, configuration_.isNullEqualNull);
    auto schema = relation->getSchema();

    for (auto col : schema->getColumns()) {
        LOG(DEBUG) << boost::format{"PLI for %1%: %2%"}
            % col->toString() % relation->getColumnData(col->getIndex())->getPositionListIndex()->toString();
    }

    auto profilingContext = std::make_shared<ProfilingContext>(
            configuration_,
            relation,
            uccConsumer_,
            fdConsumer_,
            cachingMethod_,
            evictionMethod_,
            cachingMethodValue
            );


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
    auto startTime = std::chrono::system_clock::now();
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
        totalTrickle += searchSpace->tricklingDown / 1000000;
    }
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);

    std::cout << "Time: " << elapsed_milliseconds.count() << " milliseconds" << std::endl;
    std::cout << "Error calculation count: " << totalErrorCalcCount << std::endl;
    std::cout << "Total ascension: " << totalAscension << std::endl;
    std::cout << "Total trickle: " << totalTrickle << std::endl;
    std::cout << "====RESULTS-FD====\r\n" << fdsToString();
    std::cout << "====RESULTS-UCC====\r\n" << uccsToString();
    return elapsed_milliseconds.count();
}

Pyro::Pyro(fs::path const &path, char separator, bool hasHeader, int seed, double maxError) :
        inputGenerator_(path, separator, hasHeader),
        cachingMethod_(CachingMethod::NOCACHING),
        evictionMethod_(CacheEvictionMethod::DEFAULT) {
    uccConsumer_ = [this](auto const& key) { this->discoveredUCCs_.push_back(key); };
    fdConsumer_ = [this](auto const& fd) { this->discoveredFDs_.push_back(fd); };
    configuration_.seed = seed;
    configuration_.maxUccError = maxError;
    configuration_.maxUccError = maxError;
}

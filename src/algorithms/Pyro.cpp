#include <chrono>

#include "core/FdG1Strategy.h"
#include "core/KeyG1Strategy.h"
#include "algorithms/Pyro.h"

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
    for (auto& searchSpace : searchSpaces_) {
        searchSpace->setContext(profilingContext);
        searchSpace->ensureInitialized();
        searchSpace->discover();
    }
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);

    std::cout << "Time: " << elapsed_milliseconds.count() << " milliseconds" << std::endl;
}

Pyro::Pyro(fs::path const &path) : inputGenerator_(path) {
    uccConsumer_ = [](auto const& key) { std::cout << "Found ucc: " << key.toString() << std::endl; };
    fdConsumer_ = [](auto const& key) { std::cout << "Found fd: " << key.toString() << std::endl; };
}

#include "DFD.h"

#include <boost/asio.hpp>

#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "PositionListIndex.h"
#include "LatticeTraversal/LatticeTraversal.h"

unsigned long long DFD::execute() {
    RelationalSchema const* const schema = relation->getSchema();

    if (relation->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }

    auto startTime = std::chrono::system_clock::now();

    //search for unique columns
    for (auto const& column : schema->getColumns()) {
        ColumnData& columnData = relation->getColumnData(column->getIndex());
        PositionListIndex const* const columnPLI = columnData.getPositionListIndex();

        if (columnPLI->getNumNonSingletonCluster() == 0) {
            Vertical const lhs = Vertical(*column);
            uniqueColumns.push_back(lhs);
            //we do not register an FD at once, because we check for FDs with empty LHS later
        }
    }

    double progressStep = 100.0 / schema->getNumColumns();
    boost::asio::thread_pool searchSpacePool(numberOfThreads);

    for (auto & rhs : schema->getColumns()) {
        boost::asio::post(searchSpacePool, [this, &rhs, schema, &progressStep]() {
            ColumnData const& rhsData = relation->getColumnData(rhs->getIndex());
            PositionListIndex const *const rhsPLI = rhsData.getPositionListIndex();

            /* if all the rows have the same value, then we register FD with empty LHS
             * if we have minimal FD like []->RHS, it is impossible to find smaller FD with this RHS,
             * so we register it and move to the next RHS
             * */
            if (rhsPLI->getNepAsLong() == relation->getNumTuplePairs()) {
                registerFD(*(schema->emptyVertical), *rhs);
                addProgress(progressStep);
                return;
            }

            auto searchSpace = LatticeTraversal(rhs.get(), relation.get(), uniqueColumns, partitionStorage.get());
            auto const minimalDeps = searchSpace.findLHSs();

            for (auto const& minimalDependencyLHS: minimalDeps) {
                registerFD(minimalDependencyLHS, *rhs);
            }
            addProgress(progressStep);
            std::cout << static_cast<int>(getProgress().second) << "%\n";
        });
    }

    searchSpacePool.join();
    setProgress(100);

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    long long aprioriMillis = elapsed_milliseconds.count();

    //std::cout << "====JSON-FD========\r\n" << getJsonFDs() << std::endl;
    std::cout << "> FD COUNT: " << fdCollection_.size() << std::endl;
    std::cout << "> HASH: " << FDAlgorithm::fletcher16() << std::endl;

    return aprioriMillis;
}

DFD::DFD(const std::filesystem::path &path, char separator, bool hasHeader, unsigned int parallelism)
        : FDAlgorithm(path, separator, hasHeader),
          numberOfThreads(parallelism <= 0 ? std::thread::hardware_concurrency() : parallelism) {
    relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    partitionStorage = std::make_unique<PartitionStorage>(relation.get(), CachingMethod::ALLCACHING, CacheEvictionMethod::MEDAINUSAGE);
}

void DFD::registerFD(Vertical vertical, Column rhs) {
    std::scoped_lock lock(registerFdMutex);
    FDAlgorithm::registerFD(vertical, rhs);
};

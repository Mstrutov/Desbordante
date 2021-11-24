#include "DFD.h"

#include <boost/asio.hpp>

#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "PositionListIndex.h"
#include "LatticeTraversal/LatticeTraversal.h"

unsigned long long DFD::executeInternal() {
    partitionStorage = std::make_unique<PartitionStorage>(relation_.get(),
                                                          CachingMethod::ALLCACHING,
                                                          CacheEvictionMethod::MEDAINUSAGE);
    RelationalSchema const* const schema = relation_->getSchema();

    auto startTime = std::chrono::system_clock::now();

    //search for unique columns
    for (auto const& column : schema->getColumns()) {
        ColumnData& columnData = relation_->getColumnData(column->getIndex());
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
            ColumnData const& rhsData = relation_->getColumnData(rhs->getIndex());
            PositionListIndex const *const rhsPLI = rhsData.getPositionListIndex();

            /* if all the rows have the same value, then we register FD with empty LHS
             * if we have minimal FD like []->RHS, it is impossible to find smaller FD with this RHS,
             * so we register it and move to the next RHS
             * */
            if (rhsPLI->getNepAsLong() == relation_->getNumTuplePairs()) {
                registerFD(*(schema->emptyVertical), *rhs);
                addProgress(progressStep);
                return;
            }

            auto searchSpace = LatticeTraversal(rhs.get(), relation_.get(), uniqueColumns, partitionStorage.get());
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
          numberOfThreads(parallelism <= 0 ? std::thread::hardware_concurrency() : parallelism) {}


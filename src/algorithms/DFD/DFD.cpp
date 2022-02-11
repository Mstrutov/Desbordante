#include "DFD.h"

#include <boost/asio.hpp>
#include <easylogging++.h>

#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "PositionListIndex.h"
#include "LatticeTraversal/LatticeTraversal.h"


unsigned long long DFD::ExecuteInternal() {
    partition_storage_ = std::make_unique<PartitionStorage>(relation_.get(),
                                                            CachingMethod::kAllCaching,
                                                            CacheEvictionMethod::kMedainUsage);
    RelationalSchema const* const schema = relation_->GetSchema();

    auto start_time = std::chrono::system_clock::now();

    //search for unique columns
    for (auto const& column : schema->GetColumns()) {
        ColumnData& column_data = relation_->GetColumnData(column->GetIndex());
        util::PositionListIndex const* const column_pli = column_data.GetPositionListIndex();

        if (column_pli->GetNumNonSingletonCluster() == 0) {
            Vertical const lhs = Vertical(*column);
            unique_columns_.push_back(lhs);
            //we do not register an FD at once, because we check for FDs with empty LHS later
        }
    }

    double progress_step = 100.0 / schema->GetNumColumns();
    boost::asio::thread_pool search_space_pool(number_of_threads_);

    for (auto& rhs : schema->GetColumns()) {
        boost::asio::post(search_space_pool, [this, &rhs, schema, &progress_step]() {
            ColumnData const& rhs_data = relation_->GetColumnData(rhs->GetIndex());
            util::PositionListIndex const* const rhs_pli = rhs_data.GetPositionListIndex();

            /* if all the rows have the same value, then we register FD with empty LHS
             * if we have minimal FD like []->RHS, it is impossible to find smaller FD with this RHS,
             * so we register it and move to the next RHS
             * */
            if (rhs_pli->GetNepAsLong() == relation_->GetNumTuplePairs()) {
                RegisterFd(*(schema->empty_vertical_), *rhs);
                AddProgress(progress_step);
                return;
            }

            auto search_space = LatticeTraversal(rhs.get(), relation_.get(), unique_columns_,
                                                 partition_storage_.get());
            auto const minimal_deps = search_space.FindLHSs();

            for (auto const& minimal_dependency_lhs: minimal_deps) {
                RegisterFd(minimal_dependency_lhs, *rhs);
            }
            AddProgress(progress_step);
            LOG(INFO) << static_cast<int>(GetProgress().second);
        });
    }

    search_space_pool.join();
    SetProgress(100);

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    long long apriori_millis = elapsed_milliseconds.count();

    LOG(INFO) << "> FD COUNT: " << fd_collection_.size();
    LOG(INFO) << "> HASH: " << PliBasedFDAlgorithm::Fletcher16();

    return apriori_millis;
}

DFD::DFD(const std::filesystem::path& path, char separator, bool has_header,
         unsigned int parallelism)
    : PliBasedFDAlgorithm(path, separator, has_header),
      number_of_threads_(parallelism <= 0 ? std::thread::hardware_concurrency() : parallelism) {}

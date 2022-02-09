#include <chrono>
#include <mutex>
#include <thread>

#include "FdG1Strategy.h"
#include "KeyG1Strategy.h"
#include "Pyro.h"

std::mutex searchSpacesMutex;

unsigned long long Pyro::ExecuteInternal() {
    using std::cout;
    auto start_time = std::chrono::system_clock::now();

    auto schema = relation_->GetSchema();

    auto profiling_context = std::make_unique<ProfilingContext>(
        configuration_, relation_.get(), ucc_consumer_, fd_consumer_, caching_method_,
        eviction_method_, caching_method_value_);

    std::function<bool(DependencyCandidate const&, DependencyCandidate const&)> launch_pad_order;
    if (configuration_.launch_pad_order == "arity") {
        launch_pad_order = DependencyCandidate::FullArityErrorComparator;
    } else if (configuration_.launch_pad_order == "error") {
        launch_pad_order = DependencyCandidate::FullErrorArityComparator;
    } else {
        throw std::runtime_error("Unknown comparator type");
    }

    int next_id = 0;
    if (configuration_.is_find_keys) {
        std::unique_ptr<DependencyStrategy> strategy;
        if (configuration_.ucc_error_measure == "g1prime") {
            strategy = std::make_unique<KeyG1Strategy>(configuration_.max_ucc_error,
                                                       configuration_.error_dev);
        } else {
            throw std::runtime_error("Unknown key error measure.");
        }
        search_spaces_.push_back(std::make_unique<SearchSpace>(next_id++, std::move(strategy),
                                                               schema, launch_pad_order));
    }
    if (configuration_.is_find_fds) {
        for (auto& rhs : schema->GetColumns()) {
            std::unique_ptr<DependencyStrategy> strategy;
            if (configuration_.ucc_error_measure == "g1prime") {
                strategy = std::make_unique<FdG1Strategy>(rhs.get(), configuration_.max_ucc_error,
                                                          configuration_.error_dev);
            } else {
                throw std::runtime_error("Unknown key error measure.");
            }
            search_spaces_.push_back(std::make_unique<SearchSpace>(next_id++, std::move(strategy),
                                                                   schema, launch_pad_order));
        }
    }
    unsigned long long init_time_millis = std::chrono::duration_cast<std::chrono::milliseconds>(
                                              std::chrono::system_clock::now() - start_time)
                                              .count();

    start_time = std::chrono::system_clock::now();
    unsigned int total_error_calc_count = 0;
    unsigned long long total_ascension = 0;
    unsigned long long total_trickle = 0;
    double progress_step = 100.0 / search_spaces_.size();

    const auto work_on_search_space = [this, &progress_step](
        std::list<std::unique_ptr<SearchSpace>>& search_spaces,
        ProfilingContext* profiling_context, int id) {
        unsigned long long millis = 0;
        while (true) {
            auto thread_start_time = std::chrono::system_clock::now();
            std::unique_ptr<SearchSpace> polled_space;
            {
                std::scoped_lock<std::mutex> lock(searchSpacesMutex);
                if (search_spaces.empty()) {
                    break;
                }
                polled_space = std::move(search_spaces.front());
                search_spaces.pop_front();
            }
            LOG(TRACE) << "Thread" << id << " got SearchSpace";
            polled_space->SetContext(profiling_context);
            polled_space->EnsureInitialized();
            polled_space->Discover();
            AddProgress(progress_step);

            millis += std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::system_clock::now() - thread_start_time)
                          .count();
        }
        //cout << "Thread" << id << " stopped working, ELAPSED TIME: " << millis << "ms.\n";
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < configuration_.parallelism; i++) {
        //std::thread();
        threads.emplace_back(work_on_search_space, std::ref(search_spaces_),
                             profiling_context.get(), i);
    }

    for (int i = 0; i < configuration_.parallelism; i++) {
        threads[i].join();
    }

    SetProgress(100);
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);

    LOG(DEBUG) << boost::format{"FdG1 error calculation: %1% ms"} %
                      (FdG1Strategy::nanos_ / 1000000);
    std::cout << "Init time: " << init_time_millis << "ms" << std::endl;
    std::cout << "Time: " << elapsed_milliseconds.count() << " milliseconds" << std::endl;
    std::cout << "Error calculation count: " << total_error_calc_count << std::endl;
    std::cout << "Total ascension time: " << total_ascension << "ms" << std::endl;
    std::cout << "Total trickle time: " << total_trickle << "ms" << std::endl;
    std::cout << "Total intersection time: "
              << util::PositionListIndex::micros_ / 1000 << "ms" << std::endl;
    /*std::cout << "====RESULTS-FD====\r\n" << FDsToString();
    std::cout << "====RESULTS-UCC====\r\n" << UCCsToString();
    std::cout << "====JSON-FD========\r\n" << FDAlgorithm::GetJsonFDs() << std::endl;*/

    std::cout << "HASH: " << PliBasedFDAlgorithm::Fletcher16() << std::endl;
    return elapsed_milliseconds.count();
}

Pyro::Pyro(std::filesystem::path const& path, char separator, bool has_header, int seed,
           double max_error, unsigned int max_lhs, int parallelism)
    : PliBasedFDAlgorithm(path, separator, has_header),
      caching_method_(CachingMethod::kCoin),
      eviction_method_(CacheEvictionMethod::kDefault) {
    ucc_consumer_ = [this](auto const& key) {
        this->DiscoverUcc(key);
    };
    fd_consumer_ = [this](auto const& fd) {
        this->DiscoverFd(fd);
        this->FDAlgorithm::RegisterFd(fd.lhs_, fd.rhs_);
    };
    configuration_.seed = seed;
    configuration_.max_ucc_error = max_error;
    configuration_.max_ucc_error = max_error;
    configuration_.max_lhs = max_lhs;
    configuration_.parallelism = parallelism <= 0 ?
                                 std::thread::hardware_concurrency() :
                                 parallelism;
}

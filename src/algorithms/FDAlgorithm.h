#pragma once

#include <filesystem>
#include <list>
#include <mutex>

#include <boost/any.hpp>

#include "CSVParser.h"
#include "FD.h"
#include "primitive.h"

namespace util {
class AgreeSetFactory;
}

/* It is highly recommended to inherit your Algorithm from this class.
 * Consider TANE as an example of such a FDAlgorithm usage.
 * */
class FDAlgorithm : public algos::Primitive {
public:
    /* Algorithm configuration struct */
    struct Config {
        using ParamValue = boost::any;
        using ParamsMap = std::unordered_map<std::string, ParamValue>;

        std::filesystem::path data{};   /* Path to input file */
        char separator = ',';           /* Separator for csv */
        bool has_header = true;         /* Indicates if input file has header */
        bool is_null_equal_null = true; /* Is NULL value equals another NULL value */
        unsigned int max_lhs = -1;      /* Maximum size of lhs value in fds to mine */
        ushort parallelism = 0;         /* Number of threads to use. If 0 is specified
                                         * this value is initialized to
                                         * std::thread::hardware_concurrency.
                                         */

        ParamsMap special_params{}; /* Other special parameters unique for a particular algorithm
                                     * algorithm. Use GetSpecialParam() to retrieve parameters by
                                     * name.
                                     */
    };

private:
    friend util::AgreeSetFactory;

    std::mutex mutable register_mutex_;

    void InitConfigParallelism();

protected:
    /* Algorithm configuration */
    Config config_;
    /* содержит множество найденных функциональных зависимостей. Это поле будет использоваться при тестировании,
     * поэтому важно положить сюда все намайненные ФЗ
     * */
    std::list<FD> fd_collection_;

    virtual void Initialize() = 0;
    // Main logic of the algorithm
    virtual unsigned long long ExecuteInternal() = 0;

    template <typename ParamType>
    ParamType GetSpecialParam(std::string const& param_name) const {
        Config::ParamValue const& value = config_.special_params.at(param_name);
        return boost::any_cast<ParamType>(value);
    }

public:
    constexpr static std::string_view kDefaultPhaseName = "FD mining";

    explicit FDAlgorithm(Config const& config, std::vector<std::string_view> phase_names)
        : Primitive(config.data, config.separator, config.has_header, std::move(phase_names)),
          config_(config) {
        InitConfigParallelism();
    }

    /* эти методы кладут зависимость в хранилище - можно пользоваться ими напрямую или override-нуть,
     * если нужно какое-то кастомное поведение
     * */
    virtual void RegisterFd(Vertical lhs, Column rhs) {
        std::scoped_lock lock(register_mutex_);
        fd_collection_.emplace_back(std::move(lhs), std::move(rhs));
    }
    virtual void RegisterFd(FD fd_to_register) {
        std::scoped_lock lock(register_mutex_);
        fd_collection_.push_back(std::move(fd_to_register));
    }

    /* fd_collection_ getter and setter */
    std::list<FD> const& FdList() const noexcept {
        return fd_collection_;
    }
    std::list<FD>& FdList() noexcept {
        return fd_collection_;
    }

    /* возвращает набор ФЗ в виде JSON-а. По сути, это просто представление фиксированного формата для сравнения
     * результатов разных алгоритмов. JSON - на всякий случай, если потом, например, понадобится загрузить список в
     * питон и как-нибудь его поанализировать
     * */
    std::string GetJsonFDs();

    /* Returns a vector of columns containing only unique values (i.e. keys).
     * Should be called after execute() only.
     * NOTE: retrieves keys from mined fds, so could be quite slow on wide
     * tables with many fds.
     * If your algorithm is inherited from FDAlgorithm but not from
     * PliBasedFDAlgorithm and generates ColumnLayoutRelationData from the
     * input table or in some similar way parses table, override this method
     * and use parsed table representation to retrieve keys (for performance
     * purposes).
     * PliBasedFDAlgorithm::GetKeys() is already overriding this method.
     */
    virtual std::vector<Column const*> GetKeys() const;

    // считает контрольную сумму Флетчера - нужно для тестирования по хешу
    unsigned int Fletcher16();

    unsigned long long Execute();

    virtual ~FDAlgorithm() = default;
};

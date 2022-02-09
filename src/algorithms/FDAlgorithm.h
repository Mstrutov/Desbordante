#pragma once

#include <filesystem>
#include <list>
#include <mutex>

#include "CSVParser.h"
#include "FD.h"

namespace util {
class AgreeSetFactory;
}

/* It is highly recommended to inherit your Algorithm from this class.
 * Consider TANE as an example of such a FDAlgorithm usage.
 * */
class FDAlgorithm {
private:
    friend util::AgreeSetFactory;

    std::mutex mutable register_mutex_;

    std::mutex mutable progress_mutex_;
    double cur_phase_progress_ = 0;
    uint8_t cur_phase_id_ = 0;

protected:

    /* создаётся в конструкторе, дальше предполагается передать его один раз в
     * ColumnLayoutRelationData в execute(), и больше не трогать
     * */
    CSVParser input_generator_;
    /* содержит множество найденных функциональных зависимостей. Это поле будет использоваться при тестировании,
     * поэтому важно положить сюда все намайненные ФЗ
     * */
    std::list<FD> fd_collection_;
    /* Vector of names of algorithm phases, should be initialized in a constructor
     * if algorithm has more than one phase. This vector is used to determine the
     * total number of phases
     */
    std::vector<std::string_view> const phase_names_;

    bool const is_null_equal_null_;

    void AddProgress(double const val) noexcept;
    void SetProgress(double const val) noexcept;
    void ToNextProgressPhase() noexcept;

    virtual void Initialize() = 0;
    // Main logic of the algorithm
    virtual unsigned long long ExecuteInternal() = 0;
public:
    constexpr static double kTotalProgressPercent = 100.0;

    explicit FDAlgorithm(std::filesystem::path const& path, char separator = ',',
                         bool has_header = true, bool const is_null_equal_null = true,
                         std::vector<std::string_view> phase_names = {"FD mining"})
        : input_generator_(path, separator, has_header),
          phase_names_(std::move(phase_names)),
          is_null_equal_null_(is_null_equal_null) {}

    /* эти методы кладут зависимость в хранилище - можно пользоваться ими напрямую или override-нуть,
     * если нужно какое-то кастомное поведение
     * */
    virtual void RegisterFd(Vertical lhs, Column rhs) {
        std::scoped_lock lock(progress_mutex_);
        fd_collection_.emplace_back(std::move(lhs), std::move(rhs));
    }
    virtual void RegisterFd(FD fd_to_register) {
        std::scoped_lock lock(progress_mutex_);
        fd_collection_.push_back(std::move(fd_to_register));
    }

    // геттер к набору ФЗ - нужно для тестирования
    std::list<FD> const& FdList() const { return fd_collection_; }

    /* возвращает набор ФЗ в виде JSON-а. По сути, это просто представление фиксированного формата для сравнения
     * результатов разных алгоритмов. JSON - на всякий случай, если потом, например, понадобится загрузить список в
     * питон и как-нибудь его поанализировать
     * */
    std::string GetJsonFDs();
    std::vector<std::string_view> const& GetPhaseNames() const noexcept {
        return phase_names_;
    }
    /* Returns pair with current progress state.
     * Pair has the form <current phase id, current phase progess>
     */
    std::pair<uint8_t, double> GetProgress() const noexcept;

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

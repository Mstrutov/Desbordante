#pragma once

#include <filesystem>
#include <list>
#include <mutex>

#include "CSVParser.h"
#include "FD.h"
#include "ColumnLayoutRelationData.h"

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
    uint8_t cur_phase_id = 0;

    /* создаётся в конструкторе, дальше предполагается передать его один раз в
     * ColumnLayoutRelationData в execute(), и больше не трогать
     * */
    CSVParser inputGenerator_;
protected:
    /* содержит множество найденных функциональных зависимостей. Это поле будет использоваться при тестировании,
     * поэтому важно положить сюда все намайненные ФЗ
     * */
    std::list<FD> fdCollection_;
    /* Vector of names of algorithm phases, should be initialized in a constructor
     * if algorithm has more than one phase. This vector is used to determine the
     * total number of phases
     */
    std::vector<std::string_view> const phase_names_;

    std::unique_ptr<ColumnLayoutRelationData> relation_;
    bool const is_null_equal_null_;

    void addProgress(double const val) noexcept;
    void setProgress(double const val) noexcept;
    void toNextProgressPhase() noexcept;

    // Main logic of the algorithm
    virtual unsigned long long executeInternal() = 0;
public:
    constexpr static double kTotalProgressPercent = 100.0;

    explicit FDAlgorithm (std::filesystem::path const& path,
                          char separator = ',', bool hasHeader = true,
                          bool const is_null_equal_null = true,
                          std::vector<std::string_view> phase_names = { "FD mining" })
            : inputGenerator_(path, separator, hasHeader),
              phase_names_(std::move(phase_names)), is_null_equal_null_(is_null_equal_null) {}

    /* эти методы кладут зависимость в хранилище - можно пользоваться ими напрямую или override-нуть,
     * если нужно какое-то кастомное поведение
     * */
    virtual void registerFD(Vertical lhs, Column rhs) {
        std::scoped_lock lock(progress_mutex_);
        fdCollection_.emplace_back(std::move(lhs), std::move(rhs));
    }
    virtual void registerFD(FD fdToRegister) {
        std::scoped_lock lock(progress_mutex_);
        fdCollection_.push_back(std::move(fdToRegister));
    }

    // геттер к набору ФЗ - нужно для тестирования
    std::list<FD> const& fdList() const { return fdCollection_; }

    /* возвращает набор ФЗ в виде JSON-а. По сути, это просто представление фиксированного формата для сравнения
     * результатов разных алгоритмов. JSON - на всякий случай, если потом, например, понадобится загрузить список в
     * питон и как-нибудь его поанализировать
     * */
    std::string getJsonFDs();
    std::vector<std::string_view> const& getPhaseNames() const noexcept {
        return phase_names_;
    }
    /* Returns pair with current progress state.
     * Pair has the form <current phase id, current phase progess>
     */
    std::pair<uint8_t, double> getProgress() const noexcept;

    // считает контрольную сумму Флетчера - нужно для тестирования по хешу
    unsigned int fletcher16();

    unsigned long long execute();

    ColumnLayoutRelationData const& getRelation() const noexcept {
        // getRealtion should be called after input file is parsed i.e. after algorithm execution
        assert(relation_ != nullptr);
        return *relation_;
    }

    virtual ~FDAlgorithm() = default;
};

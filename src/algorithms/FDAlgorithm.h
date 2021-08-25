#pragma once

#include <filesystem>
#include <list>
#include <atomic>

#ifndef __cpp_lib_atomic_float
#include <mutex>
#endif

#include "CSVParser.h"
#include "FD.h"

class AgreeSetFactory;

/* It is highly recommended to inherit your Algorithm from this class.
 * Consider TANE as an example of such a FDAlgorithm usage.
 * */
class FDAlgorithm {
private:
    friend AgreeSetFactory;

#ifdef __cpp_lib_atomic_float
    std::atomic<double> progress_ = 0;
#else
    double progress_ = 0;
    std::mutex mutable progress_mutex_;
#endif
protected:
    /* создаётся в конструкторе, дальше предполагается передать его один раз в
     * ColumnLayoutRelationData в execute(), и больше не трогать
     * */
    CSVParser inputGenerator_;

    /* содержит множество найденных функциональных зависимостей. Это поле будет использоваться при тестировании,
     * поэтому важно положить сюда все намайненные ФЗ
     * */
    std::list<FD> fdCollection_;

    void addProgress(double const val) noexcept;
    void setProgress(double const val) noexcept;
public:
    explicit FDAlgorithm (std::filesystem::path const& path,
                          char separator = ',', bool hasHeader = true)
            : inputGenerator_(path, separator, hasHeader) {}

    /* эти методы кладут зависимость в хранилище - можно пользоваться ими напрямую или override-нуть,
     * если нужно какое-то кастомное поведение
     * */
    virtual void registerFD(Vertical lhs, Column rhs) {
        fdCollection_.emplace_back(std::move(lhs), std::move(rhs));
    }
    virtual void registerFD(FD fdToRegister) {
        fdCollection_.push_back(std::move(fdToRegister));
    }

    // геттер к набору ФЗ - нужно для тестирования
    std::list<FD> const& fdList() const { return fdCollection_; }

    /* возвращает набор ФЗ в виде JSON-а. По сути, это просто представление фиксированного формата для сравнения
     * результатов разных алгоритмов. JSON - на всякий случай, если потом, например, понадобится загрузить список в
     * питон и как-нибудь его поанализировать
     * */
    std::string getJsonFDs();
    double getProgress() const noexcept;

    // считает контрольную сумму Флетчера - нужно для тестирования по хешу
    unsigned int fletcher16();

    // здесь будет основная логика - по сути, единстенный метод, который необходимо реализовать
    virtual unsigned long long execute() = 0;

    virtual ~FDAlgorithm() = default;
};

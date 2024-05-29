#pragma once

#include <filesystem>
#include <list>
#include <mutex>

#include <boost/any.hpp>

#include "algorithms/algorithm.h"
#include "algorithms/fd/fd.h"
#include "config/max_lhs/type.h"
#include "util/primitive_collection.h"

namespace model {
class AgreeSetFactory;
}

namespace algos {

/* It is highly recommended to inherit your Algorithm from this class.
 * Consider TANE as an example of such a FDAlgorithm usage.
 * */
class FDAlgorithm : public Algorithm {
private:
    friend model::AgreeSetFactory;

    void RegisterOptions();

    void ResetState() final;
    virtual void MakeExecuteOptsAvailableFDInternal(){};
    void MakeExecuteOptsAvailable() override;
    virtual void ResetStateFd() = 0;

protected:
    config::MaxLhsType max_lhs_;

    /* Collection of all discovered FDs
     * Every FD mining algorithm should place discovered dependecies here. Don't add new FDs by
     * accessing this field directly, use RegisterFd methods instead
     */
    util::PrimitiveCollection<FD> fd_collection_;

    /* Registers new FD.
     * Should be overrided if custom behavior is needed
     */
    virtual void RegisterFd(Vertical lhs, Column rhs) {
        if (lhs.GetArity() <= max_lhs_) fd_collection_.Register(std::move(lhs), std::move(rhs));
    }

    virtual void RegisterFd(FD fd_to_register) {
        if (fd_to_register.GetLhs().GetArity() <= max_lhs_)
            fd_collection_.Register(std::move(fd_to_register));
    }

public:
    constexpr static std::string_view kDefaultPhaseName = "FD mining";

    explicit FDAlgorithm(std::vector<std::string_view> phase_names);

    /* Returns the list of discovered FDs */
    std::list<FD> const& FdList() const noexcept {
        return fd_collection_.AsList();
    }

    std::list<FD>& FdList() noexcept {
        return fd_collection_.AsList();
    }

    std::list<FD>& SortedFdList();

    /* возвращает набор ФЗ в виде JSON-а. По сути, это просто представление фиксированного формата
     * для сравнения результатов разных алгоритмов. JSON - на всякий случай, если потом, например,
     * понадобится загрузить список в питон и как-нибудь его поанализировать
     * */
    std::string GetJsonFDs() const;

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

    virtual ~FDAlgorithm() = default;

    template <typename Container>
    static std::string FDsToJson(Container const& fds) {
        std::string result = "{\"fds\": [";
        std::vector<std::string> discovered_fd_strings;
        for (FD const& fd : fds) {
            discovered_fd_strings.push_back(fd.ToJSONString());
        }
        std::sort(discovered_fd_strings.begin(), discovered_fd_strings.end());
        for (std::string const& fd : discovered_fd_strings) {
            result += fd + ",";
        }
        if (result.back() == ',') {
            result.erase(result.size() - 1);
        }
        result += "]}";
        return result;
    }
};

}  // namespace algos

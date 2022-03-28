#pragma once

#include "PliBasedFDAlgorithm.h"
#include "custom/CustomHashes.h"

namespace algos {

class FunQuadruple {
private:
    Vertical candidate_;
    unsigned long count_;
    Vertical quasiclosure_;
    Vertical closure_;

public:
    explicit FunQuadruple(Vertical const& candidate)
        : candidate_(candidate),
          count_(0),
          quasiclosure_(*candidate.GetSchema()->empty_vertical_),
          closure_(*candidate.GetSchema()->empty_vertical_) {}
    explicit FunQuadruple(Column const& candidate) : FunQuadruple(Vertical(candidate)) {}

    Vertical const& GetCandidate() const {
        return candidate_;
    }

    unsigned long GetCount() const {
        return count_;
    }

    Vertical const& GetClosure() const {
        return closure_;
    }

    Vertical const& GetQuasiclosure() const {
        return quasiclosure_;
    }

    void SetCount(unsigned long new_count) {
        count_ = new_count;
    }

    void SetClosure(Vertical const& new_closure) {
        closure_ = new_closure;
    }

    void SetQuasiclosure(Vertical const& new_quasiclosure) {
        quasiclosure_ = new_quasiclosure;
    }

    bool operator==(FunQuadruple const& that) const {
        return candidate_ == that.candidate_;
    }

    bool operator!=(FunQuadruple const& that) const {
        return candidate_ != that.candidate_;
    }

    bool operator<(FunQuadruple const& that) const {
        return candidate_ < that.candidate_;
    }

    bool operator>(FunQuadruple const& that) const {
        return candidate_ > that.candidate_;
    }

    FunQuadruple Union(Column const& that) const;

    FunQuadruple Union(Vertical const& that) const;

    bool Contains(FunQuadruple const& that) const;

    bool Contains(Vertical const& that) const;
};

class FUN : public PliBasedFDAlgorithm {
public:
    explicit FUN(Config const& config) : PliBasedFDAlgorithm(config, {kDefaultPhaseName}) {}
    explicit FUN(std::shared_ptr<ColumnLayoutRelationData> relation, Config const& config)
        : PliBasedFDAlgorithm(std::move(relation), config, {kDefaultPhaseName}) {}

    // Entities from the algorithm itself
private:
    Vertical r_;
    Vertical r_prime_;

    using Level = std::list<FunQuadruple>;

    unsigned long long ExecuteInternal() override;

    Level GenerateCandidate(Level const& l_k) const;

    void ComputeClosure(Level& l_k_minus_1, Level const& l_k) const;

    unsigned long Count(Vertical const& l) const;

    unsigned long FastCount(Level const& l_k_minus_1, Level const& l_k,
                            FunQuadruple const& l) const;

    void ComputeQuasiClosure(Level const& l_k_minus_1, Level& l_k) const;

    void PurePrune(Level const& l_k_minus_1, Level& l_k) const;

    void DisplayFD(Level const& l_k_minus_1);

    // Supporting entities
private:
    RelationalSchema const* schema_;
    std::unordered_map<Column, std::set<Vertical>> fds_;

    bool IsKey(FunQuadruple const& l) const;
};

}  // namespace algos

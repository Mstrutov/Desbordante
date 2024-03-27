#include "canonical_od.h"

#include <sstream>

namespace algos::fastod {

template <bool Ascending>
CanonicalOD<Ascending>::CanonicalOD(AttributeSet const& context, model::ColumnIndex left,
                                    model::ColumnIndex right)
    : context_(std::move(context)), ap_(left, right) {}

template <bool Ascending>
bool CanonicalOD<Ascending>::IsValid(std::shared_ptr<DataFrame> data, PartitionCache& cache) const {
    return !(cache.GetStrippedPartition(context_, data).Swap<Ascending>(ap_.left, ap_.right));
}

template <bool Ascending>
std::string CanonicalOD<Ascending>::ToString() const {
    std::stringstream result;

    result << context_.ToString() << " : " << ap_.left + 1 << (Ascending ? "<=" : ">=") << " ~ "
           << ap_.right + 1 << "<=";

    return result.str();
}

SimpleCanonicalOD::SimpleCanonicalOD() : right_(0) {}

SimpleCanonicalOD::SimpleCanonicalOD(AttributeSet const& context, model::ColumnIndex right)
    : context_(context), right_(right) {}

bool SimpleCanonicalOD::IsValid(std::shared_ptr<DataFrame> data, PartitionCache& cache) const {
    return !(cache.GetStrippedPartition(context_, data).Split(right_));
}

std::string SimpleCanonicalOD::ToString() const {
    std::stringstream result;
    result << context_.ToString() << " : [] -> " << right_ + 1 << "<=";

    return result.str();
}

bool operator==(CanonicalOD<true> const& x, CanonicalOD<true> const& y) {
    return x.context_ == y.context_ && x.ap_ == y.ap_;
}

bool operator!=(CanonicalOD<true> const& x, CanonicalOD<true> const& y) {
    return !(x == y);
}

bool operator<(CanonicalOD<true> const& x, CanonicalOD<true> const& y) {
    if (x.ap_ != y.ap_) {
        return x.ap_ < y.ap_;
    }

    return x.context_ < y.context_;
}

bool operator==(CanonicalOD<false> const& x, CanonicalOD<false> const& y) {
    return x.context_ == y.context_ && x.ap_ == y.ap_;
}

bool operator!=(CanonicalOD<false> const& x, CanonicalOD<false> const& y) {
    return !(x == y);
}

bool operator<(CanonicalOD<false> const& x, CanonicalOD<false> const& y) {
    if (x.ap_ != y.ap_) {
        return x.ap_ < y.ap_;
    }

    return x.context_ < y.context_;
}

bool operator==(SimpleCanonicalOD const& x, SimpleCanonicalOD const& y) {
    return x.context_ == y.context_ && x.right_ == y.right_;
}

bool operator!=(SimpleCanonicalOD const& x, SimpleCanonicalOD const& y) {
    return !(x == y);
}

bool operator<(SimpleCanonicalOD const& x, SimpleCanonicalOD const& y) {
    if (x.right_ != y.right_) {
        return x.right_ < y.right_;
    }

    return x.context_ < y.context_;
}

template class CanonicalOD<true>;
template class CanonicalOD<false>;

}  // namespace algos::fastod

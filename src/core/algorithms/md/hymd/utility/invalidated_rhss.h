#pragma once

#include <iterator>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/rhss.h"
#include "model/index.h"

namespace algos::hymd::utility {
class InvalidatedRhss {
    using NewBounds = std::vector<model::md::DecisionBoundary>;

    Rhss invalidated_;
    NewBounds new_bounds_;

public:
    class UpdateView {
        InvalidatedRhss const& invalidated_info_;

        UpdateView(InvalidatedRhss const& invalidated_info) noexcept
            : invalidated_info_(invalidated_info) {}

    public:
        class update_iterator {
            using NewBoundsIterator = NewBounds::const_iterator;
            using InvalidatedIterator = Rhss::const_iterator;

            InvalidatedIterator invalidated_iter_;
            NewBoundsIterator new_bounds_iter_;

            update_iterator(InvalidatedIterator invalidated_iter,
                            NewBoundsIterator new_bounds_iter) noexcept
                : invalidated_iter_(invalidated_iter), new_bounds_iter_(new_bounds_iter) {}

        public:
            using iterator_category = std::forward_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = MdElement;
            using pointer = void;
            using reference = void;

            update_iterator operator++() noexcept {
                ++invalidated_iter_;
                ++new_bounds_iter_;
                return *this;
            }

            update_iterator operator++(int) noexcept {
                update_iterator old = *this;
                ++(*this);
                return old;
            }

            friend bool operator==(update_iterator const& iter1, update_iterator const& iter2) {
                return iter1.invalidated_iter_ == iter2.invalidated_iter_;
            }

            friend bool operator!=(update_iterator const& iter1, update_iterator const& iter2) {
                return !(iter1 == iter2);
            }

            value_type operator*() const noexcept {
                return {invalidated_iter_->index, *new_bounds_iter_};
            }

            friend UpdateView;
        };

        update_iterator begin() const noexcept {
            return {invalidated_info_.invalidated_.begin(), invalidated_info_.new_bounds_.begin()};
        }

        update_iterator end() const noexcept {
            return {invalidated_info_.invalidated_.end(), invalidated_info_.new_bounds_.end()};
        }

        friend InvalidatedRhss;
    };

    void PushBack(MdElement old_rhs, model::md::DecisionBoundary new_bound) {
        invalidated_.push_back(old_rhs);
        new_bounds_.push_back(new_bound);
    }

    UpdateView GetUpdateView() const noexcept {
        return {*this};
    }

    Rhss const& GetInvalidated() const noexcept {
        return invalidated_;
    }

    bool IsEmpty() const noexcept {
        return invalidated_.empty();
    }

    std::size_t Size() const noexcept {
        return invalidated_.size();
    }
};
}  // namespace algos::hymd::utility

#pragma once

#include "RelationalSchema.h"
#include "Vertical.h"

namespace std {
    template<>
    struct hash<Vertical> {
        size_t operator()(Vertical const& k) const {
            return k.getColumnIndices().to_ulong();
        }
    };

    template<>
    struct hash<Column> {
        size_t operator()(Column const& k) const {
            boost::dynamic_bitset<> columnIndex(k.getSchema()->getNumColumns());
            columnIndex.set(k.getIndex());
            return columnIndex.to_ulong();
        }
    };

    template<>
    struct hash<std::shared_ptr<Vertical>> {
        size_t operator()(std::shared_ptr<Vertical> const& k) const {
            return k->getColumnIndices().to_ulong();
        }
    };

    template<>
    struct hash<std::shared_ptr<Column>> {
        size_t operator()(std::shared_ptr<Column> const& k) const {
            boost::dynamic_bitset<> columnIndex(k->getSchema()->getNumColumns());
            columnIndex.set(k->getIndex());
            return columnIndex.to_ulong();
        }
    };

    template<class T>
    struct hash<std::pair<Vertical, T>> {
        size_t operator()(std::pair<Vertical, T> const& k) const {
            return k.first.getColumnIndices().to_ulong();
        }
    };

    template<class T>
    struct hash<std::pair<std::shared_ptr<Vertical>, T>> {
        size_t operator()(std::pair<std::shared_ptr<Vertical>, T> const& k) const {
            return k.first->getColumnIndices().to_ulong();
        }
    };
}


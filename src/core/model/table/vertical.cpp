#include "vertical.h"

#include <utility>

#include "model/table/relational_schema.h"

Vertical::Vertical(RelationalSchema const* rel_schema, boost::dynamic_bitset<> indices)
    : column_indices_(std::move(indices)), schema_(rel_schema) {}

Vertical::Vertical(Column const& col) : schema_(col.GetSchema()) {
    column_indices_ = boost::dynamic_bitset<>(schema_->GetNumColumns());
    column_indices_.set(col.GetIndex());
}

bool Vertical::Contains(Vertical const& that) const {
    boost::dynamic_bitset<> const& that_indices = that.column_indices_;
    if (column_indices_.size() < that_indices.size()) return false;

    return that.column_indices_.is_subset_of(column_indices_);
}

bool Vertical::Contains(Column const& that) const {
    return column_indices_.test(that.GetIndex());
}

bool Vertical::Intersects(Vertical const& that) const {
    boost::dynamic_bitset<> const& that_indices = that.column_indices_;
    return column_indices_.intersects(that_indices);
}

Vertical Vertical::Union(Vertical const& that) const {
    boost::dynamic_bitset<> retained_column_indices(column_indices_);
    retained_column_indices |= that.column_indices_;
    return schema_->GetVertical(retained_column_indices);
}

Vertical Vertical::Union(Column const& that) const {
    boost::dynamic_bitset<> retained_column_indices(column_indices_);
    retained_column_indices.set(that.GetIndex());
    return schema_->GetVertical(retained_column_indices);
}

Vertical Vertical::Project(Vertical const& that) const {
    boost::dynamic_bitset<> retained_column_indices(column_indices_);
    retained_column_indices &= that.column_indices_;
    return schema_->GetVertical(retained_column_indices);
}

Vertical Vertical::Without(Vertical const& that) const {
    boost::dynamic_bitset<> retained_column_indices(column_indices_);
    retained_column_indices &= ~that.column_indices_;
    return schema_->GetVertical(retained_column_indices);
}

Vertical Vertical::Without(Column const& that) const {
    boost::dynamic_bitset<> retained_column_indices(column_indices_);
    retained_column_indices.reset(that.GetIndex());
    return schema_->GetVertical(retained_column_indices);
}

Vertical Vertical::Invert() const {
    boost::dynamic_bitset<> flipped_indices(column_indices_);
    flipped_indices.resize(schema_->GetNumColumns());
    flipped_indices.flip();
    return schema_->GetVertical(flipped_indices);
}

Vertical Vertical::Invert(Vertical const& scope) const {
    boost::dynamic_bitset<> flipped_indices(column_indices_);
    flipped_indices ^= scope.column_indices_;
    return schema_->GetVertical(flipped_indices);
}

std::unique_ptr<Vertical> Vertical::EmptyVertical(RelationalSchema const* rel_schema) {
    return std::make_unique<Vertical>(rel_schema,
                                      boost::dynamic_bitset<>(rel_schema->GetNumColumns()));
}

std::vector<Column const*> Vertical::GetColumns() const {
    std::vector<Column const*> columns;
    for (size_t index = column_indices_.find_first(); index != boost::dynamic_bitset<>::npos;
         index = column_indices_.find_next(index)) {
        columns.push_back(&schema_->GetColumns()[index]);
    }
    return columns;
}

std::vector<unsigned> Vertical::GetColumnIndicesAsVector() const {
    std::vector<unsigned> columns;
    for (size_t index = column_indices_.find_first(); index != boost::dynamic_bitset<>::npos;
         index = column_indices_.find_next(index)) {
        columns.push_back(schema_->GetColumns()[index].GetIndex());
    }
    return columns;
}

std::string Vertical::ToString() const {
    std::string result = "[";

    if (column_indices_.find_first() == boost::dynamic_bitset<>::npos) return "[]";

    for (size_t index = column_indices_.find_first(); index != boost::dynamic_bitset<>::npos;
         index = column_indices_.find_next(index)) {
        result += schema_->GetColumn(index)->GetName();
        if (column_indices_.find_next(index) != boost::dynamic_bitset<>::npos) {
            result += ' ';
        }
    }

    result += ']';

    return result;
}

std::string Vertical::ToIndicesString() const {
    std::string result = "[";

    if (column_indices_.find_first() == boost::dynamic_bitset<>::npos) {
        return "[]";
    }

    for (size_t index = column_indices_.find_first(); index != boost::dynamic_bitset<>::npos;
         index = column_indices_.find_next(index)) {
        result += std::to_string(index);
        if (column_indices_.find_next(index) != boost::dynamic_bitset<>::npos) {
            result += ',';
        }
    }

    result += ']';

    return result;
}

std::vector<Vertical> Vertical::GetParents() const {
    if (GetArity() < 2) return std::vector<Vertical>();
    std::vector<Vertical> parents(GetArity());
    int i = 0;
    for (size_t column_index = column_indices_.find_first();
         column_index != boost::dynamic_bitset<>::npos;
         column_index = column_indices_.find_next(column_index)) {
        auto parent_column_indices = column_indices_;
        parent_column_indices.reset(column_index);
        parents[i++] = GetSchema()->GetVertical(std::move(parent_column_indices));
    }
    return parents;
}

bool Vertical::operator<(Vertical const& rhs) const {
    assert(*schema_ == *rhs.schema_);
    if (this->column_indices_ == rhs.column_indices_) return false;

    boost::dynamic_bitset<> const& lr_xor = (this->column_indices_ ^ rhs.column_indices_);
    return rhs.column_indices_.test(lr_xor.find_first());
}

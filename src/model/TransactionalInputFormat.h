#pragma once

class InputFormat {
public:
    virtual unsigned tid_column_index() const noexcept = 0;
    virtual unsigned item_column_index() const noexcept = 0;
    virtual bool tid_presence() const noexcept = 0;

    virtual ~InputFormat() = default;
};

class Singular : public InputFormat {
private:
    unsigned tid_column_index_;
    unsigned item_column_index_;

public:
    explicit Singular(unsigned tid_column_index, unsigned item_column_index)
        :tid_column_index_(tid_column_index), item_column_index_(item_column_index) {}

    unsigned  tid_column_index() const noexcept override { return tid_column_index_; }
    unsigned item_column_index() const noexcept override { return  item_column_index_; }
    bool tid_presence() const noexcept override { return true; }
};

class Tabular : public InputFormat {
private:
    bool has_tid_;

public:
    explicit Tabular(bool has_tid)
        : has_tid_(has_tid) {}

    unsigned  tid_column_index() const noexcept override { return 0; }
    unsigned item_column_index() const noexcept override { return 0; }
    bool tid_presence() const noexcept override { return has_tid_; }
};

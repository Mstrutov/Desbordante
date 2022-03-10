#include <filesystem>
#include <functional>
#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Types.h"

namespace tests {

namespace mo = model;

template <typename T>
class TestNumeric : public ::testing::Test {
public:
    using NumericType = T;
    using UnderlyingType = typename T::UnderlyingType;

private:
    using NumericTypeMethod = std::byte* (NumericType::*)(std::byte const*, std::byte const*,
                                                          std::byte*) const;

    std::unique_ptr<std::byte[]> actual_owner_;
    std::unique_ptr<std::byte[]> literal_owner_;

protected:
    std::unique_ptr<NumericType> type_;
    std::byte* actual_ptr_;
    std::byte* literal_ptr_;

    void SetUp() override {
        type_ = std::make_unique<NumericType>();
        actual_owner_.reset(type_->MakeValue(0));
        literal_owner_.reset(type_->MakeValue(0));
        actual_ptr_ = actual_owner_.get();
        literal_ptr_ = literal_owner_.get();
    }

    void SetActualAndLiteral(UnderlyingType actual, UnderlyingType literal) {
        SetActual(actual);
        SetLiteral(literal);
    }

    void SetActual(UnderlyingType actual) {
        mo::Type::GetValue<UnderlyingType>(actual_ptr_) = actual;
    }

    void SetLiteral(UnderlyingType literal) {
        mo::Type::GetValue<UnderlyingType>(literal_ptr_) = literal;
    }

    template <typename F>
    void TestBinary(UnderlyingType l, UnderlyingType r, F result_f,
                    NumericTypeMethod method_to_test) {
        SetActualAndLiteral(l, r);
        ExpectEq(result_f(l, r),
                 std::invoke(method_to_test, type_.get(), actual_ptr_, literal_ptr_, actual_ptr_),
                 GetErrorStr(l, r));
    }

    template <typename K, typename Y>
    static std::string GetErrorStr(K l, Y r) {
        return "Failed with l=" + std::to_string(l) + " and r=" + std::to_string(r);
    }
    static void ExpectEq(UnderlyingType expected, std::byte const* v, std::string error = "") {
        UnderlyingType actual = mo::Type::GetValue<decltype(expected)>(v);
        if constexpr (std::is_floating_point_v<UnderlyingType>) {
            EXPECT_DOUBLE_EQ(expected, actual) << error;
        } else {
            EXPECT_EQ(expected, actual) << error;
        }
    };
};

using NumericTypes = ::testing::Types<mo::DoubleType, mo::IntType>;
TYPED_TEST_SUITE(TestNumeric, NumericTypes);

TYPED_TEST(TestNumeric, Compare) {
    using Type = typename TypeParam::UnderlyingType;
    auto test = [this](Type l, Type r, mo::CompareResult res) {
        TestFixture::SetActualAndLiteral(l, r);
        EXPECT_EQ(res, this->type_->Compare(this->actual_ptr_, this->literal_ptr_));
    };

    test(3, 3, mo::CompareResult::kEqual);
    test(1, 2, mo::CompareResult::kLess);
    test(100, 12, mo::CompareResult::kGreater);
    test(-3, 3, mo::CompareResult::kLess);
    test(0, -100, mo::CompareResult::kGreater);
}

TYPED_TEST(TestNumeric, Negate) {
    auto test = [this](typename TypeParam::UnderlyingType v) {
        TestFixture::SetActual(v);
        TestFixture::ExpectEq(-v, this->type_->Negate(this->actual_ptr_, this->actual_ptr_));
    };

    test(0);
    test(-123.5);
    test(321.4);
}

TYPED_TEST(TestNumeric, Abs) {
    auto test = [this](typename TypeParam::UnderlyingType v) {
        TestFixture::SetActual(v);
        TestFixture::ExpectEq(std::abs(v), this->type_->Abs(this->actual_ptr_, this->actual_ptr_));
    };

    test(0);
    test(-123.5);
    test(321.4);
}

TYPED_TEST(TestNumeric, Add) {
    using Type = typename TypeParam::UnderlyingType;
    auto test = [this](Type l, Type r) {
        this->TestBinary(l, r, std::plus<Type>(), &TypeParam::Add);
    };

    test(0, 1);
    test(0, 0);
    test(1, 1);
    test(-4, 123);
    test(-123, -11);
    test(-123, 808);
}

TYPED_TEST(TestNumeric, Div) {
    using Type = typename TypeParam::UnderlyingType;
    auto test = [this](Type l, Type r) {
        this->TestBinary(l, r, std::divides<Type>(), &TypeParam::Div);
    };

    test(0, 100);
    test(22, 1);
    test(123, 321);
    test(11.4, 3.14);
    test(-102, 11);
    test(-123, 123);
    test(-21, -7);
}

TYPED_TEST(TestNumeric, Sub) {
    using Type = typename TypeParam::UnderlyingType;
    auto test = [this](Type l, Type r) {
        this->TestBinary(l, r, std::minus<Type>(), &TypeParam::Sub);
    };

    test(0, 100);
    test(22, 12);
    test(123, 321);
    test(2.72, 1.3123141);
    test(-102, 11);
    test(-123, 123);
    test(-21, -7);
}

TYPED_TEST(TestNumeric, Mul) {
    using Type = typename TypeParam::UnderlyingType;
    auto test = [this](Type l, Type r) {
        this->TestBinary(l, r, std::multiplies<Type>(), &TypeParam::Mul);
    };

    test(0, 100);
    test(100, 0);
    test(22, 12);
    test(123, 321);
    test(2.72, 1.3123141);
    test(-102, 11);
    test(-123, 123);
    test(-21, -7);
}

TYPED_TEST(TestNumeric, Pow) {
    using Type = typename TypeParam::UnderlyingType;
    auto test = [this](Type l, double r) {
        this->SetActual(l);
        this->ExpectEq(std::pow(l, r), this->type_->Power(this->actual_ptr_, r, this->actual_ptr_),
                       TestFixture::GetErrorStr(l, r));
    };

    test(0, 100);
    test(22, 12);
    test(123, 321);
    test(2.72, 1.3123141);
    test(-102, 11);
    test(-123, 123);
    test(-21, -7);
}

TYPED_TEST(TestNumeric, ValueToString) {
    auto test = [this](typename TypeParam::UnderlyingType v) {
        this->SetActual(v);
        EXPECT_EQ(std::to_string(v), this->type_->ValueToString(this->actual_ptr_));
    };

    test(0);
    test(123);
    test(3.14123123182387);
    test(-1231.123456678987654321);
}

}  // namespace tests

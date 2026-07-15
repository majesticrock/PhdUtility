#include <complex>
#include <iostream>
#include <type_traits>

// Include your header here
#include <mrock/utility/UnderlyingRealType.hpp>

int main() {
    using mrock::utility::UnderlyingRealType;
    using mrock::utility::UnderlyingRealType_t;

    {
        std::cout << "Test 1: plain real types\n";

        constexpr bool int_is_int = std::is_same_v<UnderlyingRealType_t<int>, int>;
        constexpr bool float_is_float = std::is_same_v<UnderlyingRealType_t<float>, float>;
        constexpr bool double_is_double = std::is_same_v<UnderlyingRealType_t<double>, double>;

        static_assert(int_is_int);
        static_assert(float_is_float);
        static_assert(double_is_double);

        std::cout << "Plain real type tests passed.\n\n";
    }

    {
        std::cout << "Test 2: complex types\n";

        constexpr bool complex_float_is_float =
            std::is_same_v<UnderlyingRealType_t<std::complex<float>>, float>;

        constexpr bool complex_double_is_double =
            std::is_same_v<UnderlyingRealType_t<std::complex<double>>, double>;

        constexpr bool complex_long_double_is_long_double =
            std::is_same_v<UnderlyingRealType_t<std::complex<long double>>, long double>;

        static_assert(complex_float_is_float);
        static_assert(complex_double_is_double);
        static_assert(complex_long_double_is_long_double);

        std::cout << "Complex type tests passed.\n\n";
    }

    {
        std::cout << "Test 3: const, volatile, and reference qualifiers\n";

        constexpr bool const_double_is_double =
            std::is_same_v<UnderlyingRealType_t<const double>, double>;

        constexpr bool double_ref_is_double =
            std::is_same_v<UnderlyingRealType_t<double&>, double>;

        constexpr bool const_double_ref_is_double =
            std::is_same_v<UnderlyingRealType_t<const double&>, double>;

        constexpr bool volatile_float_ref_is_float =
            std::is_same_v<UnderlyingRealType_t<volatile float&>, float>;

        static_assert(const_double_is_double);
        static_assert(double_ref_is_double);
        static_assert(const_double_ref_is_double);
        static_assert(volatile_float_ref_is_float);

        std::cout << "cv/ref qualifier tests passed.\n\n";
    }

    {
        std::cout << "Test 4: qualified complex types\n";

        using ComplexDouble = std::complex<double>;
        using ComplexFloat = std::complex<float>;

        constexpr bool const_complex_double_is_double =
            std::is_same_v<UnderlyingRealType_t<const ComplexDouble>, double>;

        constexpr bool complex_double_ref_is_double =
            std::is_same_v<UnderlyingRealType_t<ComplexDouble&>, double>;

        constexpr bool const_complex_float_ref_is_float =
            std::is_same_v<UnderlyingRealType_t<const ComplexFloat&>, float>;

        constexpr bool volatile_complex_float_ref_is_float =
            std::is_same_v<UnderlyingRealType_t<volatile ComplexFloat&>, float>;

        static_assert(const_complex_double_is_double);
        static_assert(complex_double_ref_is_double);
        static_assert(const_complex_float_ref_is_float);
        static_assert(volatile_complex_float_ref_is_float);

        std::cout << "Qualified complex type tests passed.\n\n";
    }

    {
        std::cout << "Test 5: direct struct usage\n";

        constexpr bool struct_plain =
            std::is_same_v<typename UnderlyingRealType<int>::type, int>;

        constexpr bool struct_complex =
            std::is_same_v<typename UnderlyingRealType<std::complex<double>>::type, double>;

        static_assert(struct_plain);
        static_assert(struct_complex);

        std::cout << "Direct struct usage tests passed.\n\n";
    }

    std::cout << "All UnderlyingRealType tests passed.\n";

    return 0;
}
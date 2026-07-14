#ifndef MROCK_UTILITY_INCLUDE_MROCK_UTILITY_CONSTEXPR_POWER_HPP
#define MROCK_UTILITY_INCLUDE_MROCK_UTILITY_CONSTEXPR_POWER_HPP
namespace mrock::utility {
template <int exponent, class NumberType, class ResultType = double>
constexpr ResultType constexpr_power(NumberType base) {
    if constexpr (exponent < 0) {
        return 1. / constexpr_power<-exponent, NumberType, ResultType>(base);
    } else if constexpr (exponent > 0) {
        if constexpr (exponent & 1) {
            return base * constexpr_power<exponent - 1, NumberType, ResultType>(base);
        } else {
            return constexpr_power<exponent / 2, NumberType, ResultType>(base * base);
        }
    } else {
        return ResultType{1};
    }
}
}  // namespace mrock::utility
#endif  // MROCK_UTILITY_INCLUDE_MROCK_UTILITY_CONSTEXPR_POWER_HPP

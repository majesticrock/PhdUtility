#pragma once

namespace Utility {
    template<int exponent, class NumberType, class ResultType = double>
    constexpr ResultType constexprPower(NumberType base) 
    {
        if constexpr (exponent < 0)
        {
            return 1. / constexprPower<-exponent, NumberType, ResultType>(base);
        }
        else if constexpr (exponent > 0)
        {
            if constexpr (exponent % 2 == 0)
            {
                return constexprPower<exponent - 1, NumberType, ResultType>(base * base);
            } 
            else 
            {
                return base * constexprPower<exponent - 1, NumberType, ResultType>(base);
            }
        }
        else 
        {
            return 1;
        }
    }
}
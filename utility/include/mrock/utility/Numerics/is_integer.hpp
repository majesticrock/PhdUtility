#pragma once
#include <cmath>

namespace mrock::utility::Numerics {
    template<class T, class tol_t=double>
    constexpr bool is_integer(T val, tol_t tol = 1e-12) {
        using std::abs;
        return abs(val - (T)((int)val)) < tol || abs(val - (T)((int)val + 1)) < tol;
    }
}
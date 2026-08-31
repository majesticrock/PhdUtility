#pragma once
#include <cassert>

namespace mrock::symbolic_operators {

template <int L>
struct IntMomentum {
private:
    constexpr int wrap(int x) const noexcept {
        x %= L;
        if (x < 0) x += L;
        return x;
    }

public:
    int i = 0; // stored index in [0,L)

    constexpr IntMomentum() = default;

    constexpr IntMomentum(int idx)
        : i(wrap(idx))
    {
        assert(L > 0);
        assert(L % 2 == 0);
    }

    constexpr IntMomentum& operator+=(const IntMomentum& rhs) noexcept {
        i = wrap(i + rhs.i + L / 2);
        return *this;
    }

    constexpr IntMomentum& operator-=(const IntMomentum& rhs) noexcept {
        i = wrap(i - rhs.i - L / 2);
        return *this;
    }

    constexpr IntMomentum& operator*=(int n) noexcept {
        // Work in coordinates centered at k=0.
        i = wrap((i - L/2) * n + L/2);
        return *this;
    }

    constexpr IntMomentum operator+(IntMomentum rhs) const noexcept {
        return rhs += *this;
    }

    constexpr IntMomentum operator-(IntMomentum rhs) const noexcept {
        return IntMomentum(*this) -= rhs;
    }

    constexpr IntMomentum operator*(int n) const noexcept {
        auto tmp = *this;
        tmp *= n;
        return tmp;
    }

    constexpr IntMomentum operator-() const noexcept {
        return (*this) * -1;
    }

    constexpr bool operator==(const IntMomentum& rhs) const noexcept {
        return i == rhs.i;
    }

    constexpr operator std::size_t() const noexcept {
        return static_cast<std::size_t>(i);
    }
};

template <int L>
constexpr IntMomentum<L> operator*(int n, IntMomentum<L> k) noexcept {
    return k *= n;
}

}
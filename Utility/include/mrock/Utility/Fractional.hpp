#pragma once
#include <numeric>
#include <type_traits>
#include <iostream>

namespace mrock::utility {
    template<class _int>
    struct Fractional {
        _int numerator{};
        _int denominator{1};

        constexpr Fractional() {};
        constexpr Fractional(_int integer) : numerator(integer) {};
        constexpr Fractional(_int _numerator, _int _denominator) : numerator(_numerator), denominator(_denominator) {};

        template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& numerator;
			ar& denominator;
		}

        template<class _float = double>
        constexpr _float value() const {
            return static_cast<_float>(numerator) / static_cast<_float>(denominator);
        }
        constexpr bool is_integer() const {
            return (numerator % denominator == 0);
        }
        inline void reduce_fraction() {
            if constexpr (std::is_signed_v<_int>) {
                if(denominator < 0) {
                    numerator *= -1;
                    denominator *= -1;
                }
            }
            const _int gcd = std::gcd(numerator, denominator);
            numerator /= gcd;
            denominator /= gcd;
        }

        template<typename Integer,
            std::enable_if_t<std::is_integral<Integer>::value, bool> = true>
        constexpr operator Integer() const {
            return numerator / denominator;
        }
        template<typename Floating,
            std::enable_if_t<std::is_floating_point<Floating>::value, bool> = true>
        constexpr operator Floating() const {
            return value<Floating>();
        }

        constexpr auto operator<=>(const Fractional<_int>& other) const {
            return (this->numerator * other.denominator) <=> (other.numerator * this->denominator);
        }
        constexpr bool operator==(const Fractional<_int>& other) const {
            return (this->numerator * other.denominator) == (other.numerator * this->denominator);
        };
        constexpr bool operator!=(const Fractional<_int>& other) const {
            return !(*this == other);
        };

        constexpr auto operator<=>(_int other) const {
            return this->numerator <=> (other * this->denominator);
        }
        constexpr bool operator==(_int other) const {
            return this->numerator == (other * this->denominator);
        };
        constexpr bool operator!=(_int other) const {
            return !(*this == other);
        };

        inline Fractional& operator+=(Fractional const& other) {
            if(other.denominator == this->denominator) {
                this->numerator += other.numerator;
                this->reduce_fraction();
                return *this;
            }
            this->numerator = this->numerator * other.denominator + other.numerator * this->denominator;
            this->denominator *= other.denominator;
            this->reduce_fraction();
            return *this;
        }
        inline Fractional& operator-=(Fractional const& other) {
            if(other.denominator == this->denominator) {
                this->numerator -= other.numerator;
                this->reduce_fraction();
                return *this;
            }
            this->numerator = this->numerator * other.denominator - other.numerator * this->denominator;
            this->denominator *= other.denominator;
            this->reduce_fraction();
            return *this;
        }
        inline Fractional& operator*=(Fractional const& other) {
            this->denominator *= other.denominator;
            this->numerator *= other.numerator;
            this->reduce_fraction();
            return *this;
        }
        inline Fractional& operator/=(Fractional const& other) {
            this->denominator *= other.numerator;
            this->numerator *= other.denominator;
            this->reduce_fraction();
            return *this;
        }

        constexpr Fractional& operator+=(_int other) {
            this->numerator += other * this->denominator;
            return *this;
        }
        constexpr Fractional& operator-=(_int other) {
            this->numerator -= other * this->denominator;
            return *this;
        }
        inline Fractional& operator*=(_int other) {
            this->numerator *= other;
            this->reduce_fraction();
            return *this;
        }
        inline Fractional& operator/=(_int other) {
            this->denominator *= other;
            this->reduce_fraction();
            return *this;
        }
    };

    template<class _int>
    constexpr Fractional<_int> pow(Fractional<_int> base, int exponent) {
        if(exponent == 0) return Fractional<_int>{_int(1), _int(1)};
        if(exponent < 0) return pow(Fractional<_int>{base.denominator, base.numerator}, exponent);
        if(exponent == 1) return base;
        if(exponent & 1) {
            return base * pow(base) * pow(base);
        } 
        else {
            return pow(base) * pow(base);
        }
    }

    template<class _int>
    std::ostream& operator<<(std::ostream& os, const Fractional<_int>& frac) {
        if(frac.is_integer()) {
            os << frac.numerator / frac.denominator;
        }
        else {
            if(frac.numerator < 0) {
                os << "-";
            }
            os << "(" << std::abs(frac.numerator) << "/" << frac.denominator << ")";
        }
        return os;
    }

    template<class _int>
    inline Fractional<_int> operator+(Fractional<_int> lhs, const Fractional<_int>& rhs) {
        return lhs += rhs;
    }
    template<class _int>
    inline Fractional<_int> operator-(Fractional<_int> lhs, const Fractional<_int>& rhs) {
        return lhs -= rhs;
    }
    template<class _int>
    inline Fractional<_int> operator*(Fractional<_int> lhs, const Fractional<_int>& rhs) {
        return lhs *= rhs;
    }
    template<class _int>
    inline Fractional<_int> operator/(Fractional<_int> lhs, const Fractional<_int>& rhs) {
        return lhs /= rhs;
    }

    template<class _int>
    inline Fractional<_int> operator-(Fractional<_int> rhs) {
        return rhs *= _int{-1};
    }

    template<class _int>
    constexpr Fractional<_int> operator+(Fractional<_int> lhs, _int rhs) {
        return lhs += rhs;
    }
    template<class _int>
    constexpr Fractional<_int> operator-(Fractional<_int> lhs, _int rhs) {
        return lhs -= rhs;
    }
    template<class _int>
    inline Fractional<_int> operator*(Fractional<_int> lhs, _int rhs) {
        return lhs *= rhs;
    }
    template<class _int>
    inline Fractional<_int> operator/(Fractional<_int> lhs, _int rhs) {
        return lhs /= rhs;
    }

    template<class _int>
    constexpr Fractional<_int> operator+(_int lhs, Fractional<_int> rhs) {
        return rhs += lhs;
    }
    template<class _int>
    constexpr Fractional<_int> operator-(_int lhs, Fractional<_int> rhs) {
        return rhs -= lhs;
    }
    template<class _int>
    inline Fractional<_int> operator*(_int lhs, Fractional<_int> rhs) {
        return rhs *= lhs;
    }
    template<class _int>
    inline Fractional<_int> operator/(_int lhs, Fractional<_int> rhs) {
        return rhs /= lhs;
    }

    template<class _Number, class _int>
    inline _Number operator+(const Fractional<_int>& lhs, const _Number& rhs) {
        return rhs + static_cast<_Number>(lhs);
    }
    template<class _Number, class _int>
    inline _Number operator-(const Fractional<_int>& lhs, const _Number& rhs) {
        return rhs - static_cast<_Number>(lhs);
    }
    template<class _Number, class _int>
    inline _Number operator*(const Fractional<_int>& lhs, const _Number& rhs) {
        return rhs * static_cast<_Number>(lhs);
    }
    template<class _Number, class _int>
    inline _Number operator/(const Fractional<_int>& lhs, const _Number& rhs) {
        return rhs / static_cast<_Number>(lhs);
    }
    template<class _Number, class _int>
    inline _Number operator+(const _Number& lhs, const Fractional<_int>& rhs) {
        return lhs + static_cast<_Number>(rhs);
    }
    template<class _Number, class _int>
    inline _Number operator-(const _Number& lhs, const Fractional<_int>& rhs) {
        return lhs - static_cast<_Number>(rhs);
    }
    template<class _Number, class _int>
    inline _Number operator*(const _Number& lhs, const Fractional<_int>& rhs) {
        return lhs * static_cast<_Number>(rhs);
    }
    template<class _Number, class _int>
    inline _Number operator/(const _Number& lhs, const Fractional<_int>& rhs) {
        return lhs / static_cast<_Number>(rhs);
    }
}
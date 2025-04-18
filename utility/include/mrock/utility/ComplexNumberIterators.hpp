#pragma once
#include <complex>
#include <iterator>

namespace mrock::utility {
    template<typename T>
    class RealPartIterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = const T*;
        using reference = const T&;

        RealPartIterator(std::complex<T> const * ptr) : ptr(ptr) {}

        template <class Vector>
        static RealPartIterator Begin(const Vector& vec) {
            return RealPartIterator(vec.data());
        }
        template <class Vector>
        static RealPartIterator End(const Vector& vec) {
            return RealPartIterator(vec.data() + vec.size());
        }

        T operator*() const { return ptr->real(); }
        T operator[](difference_type n) const {
            return (ptr + n)->real();
        }

        RealPartIterator& operator++() { ++ptr; return *this; }
        RealPartIterator operator++(int) { RealPartIterator tmp = *this; ++(*this); return tmp; }
        RealPartIterator& operator--() { --ptr; return *this; }
        RealPartIterator operator--(int) { RealPartIterator tmp = *this; --(*this); return tmp; }

        RealPartIterator operator+(difference_type n) const { return RealPartIterator(ptr + n); }
        RealPartIterator& operator+=(difference_type n) { ptr += n; return *this; }
        RealPartIterator operator-(difference_type n) const { return RealPartIterator(ptr - n); }
        RealPartIterator& operator-=(difference_type n) { ptr -= n; return *this; }

        difference_type operator-(const RealPartIterator& other) const { return ptr - other.ptr; }

        auto operator<=>(RealPartIterator const& other) const = default;

    private:
        std::complex<T> const * ptr;
    };

    template<typename T>
    class ImagPartIterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using reference = T&;

        ImagPartIterator(std::complex<T> const * ptr) : ptr(ptr) {}

        template <class Vector>
        static ImagPartIterator Begin(const Vector& vec) {
            return ImagPartIterator(vec.data());
        }
        template <class Vector>
        static ImagPartIterator End(const Vector& vec) {
            return ImagPartIterator(vec.data() + vec.size());
        }

        T operator*() const { return ptr->imag(); }
        T operator[](difference_type n) const {
            return (ptr + n)->imag();
        }

        ImagPartIterator& operator++() { ++ptr; return *this; }
        ImagPartIterator operator++(int) { ImagPartIterator tmp = *this; ++(*this); return tmp; }
        ImagPartIterator& operator--() { --ptr; return *this; }
        ImagPartIterator operator--(int) { ImagPartIterator tmp = *this; --(*this); return tmp; }

        ImagPartIterator operator+(difference_type n) const { return ImagPartIterator(ptr + n); }
        ImagPartIterator& operator+=(difference_type n) { ptr += n; return *this; }
        ImagPartIterator operator-(difference_type n) const { return ImagPartIterator(ptr - n); }
        ImagPartIterator& operator-=(difference_type n) { ptr -= n; return *this; }

        difference_type operator-(const ImagPartIterator& other) const { return ptr - other.ptr; }

        auto operator<=>(ImagPartIterator const& other) const = default;

    private:
        std::complex<T> const * ptr;
    };

    template<typename T>
    RealPartIterator<T> make_real_part_iterator(std::complex<T> const * ptr) {
        return RealPartIterator<T>(ptr);
    }
    template<typename T>
    RealPartIterator<T> make_real_part_iterator_end(std::complex<T> const * ptr, std::size_t length) {
        return RealPartIterator<T>(ptr + length);
    }
    template<typename T>
    ImagPartIterator<T> make_imag_part_iterator(std::complex<T> const * ptr) {
        return ImagPartIterator<T>(ptr);
    }
    template<typename T>
    ImagPartIterator<T> make_imag_part_iterator_end(std::complex<T> const * ptr, std::size_t length) {
        return ImagPartIterator<T>(ptr + length);
    }

    template<typename Vector>
    RealPartIterator<typename Vector::value_type::value_type> make_real_part_iterator(Vector const& vec) {
        return RealPartIterator<typename Vector::value_type::value_type>::Begin(vec);
    }
    template<typename Vector>
    RealPartIterator<typename Vector::value_type::value_type> make_real_part_iterator_end(Vector const& vec) {
        return RealPartIterator<typename Vector::value_type::value_type>::End(vec);
    }
    template<typename Vector>
    ImagPartIterator<typename Vector::value_type::value_type> make_imag_part_iterator(Vector const& vec) {
        return ImagPartIterator<typename Vector::value_type::value_type>::Begin(vec);
    }
    template<typename Vector>
    ImagPartIterator<typename Vector::value_type::value_type> make_imag_part_iterator_end(Vector const& vec) {
        return ImagPartIterator<typename Vector::value_type::value_type>::End(vec);
    }
}
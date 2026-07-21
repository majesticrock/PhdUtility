#include <boost/math/tools/minima.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

constexpr auto pyarray_flags = py::array::c_style | py::array::forcecast;

using py_array_double = py::array_t<double, pyarray_flags>;
using py_array_cmplx = py::array_t<std::complex<double>, pyarray_flags>;

namespace detail {

constexpr int MAX_EFFECTIVE_ROOT_TOL_BITS = 26;
constexpr double BOOST_SHIFT = 1.0;
constexpr double EPS = 1e-6;
constexpr double MIN_SEARCH_RANGE = 1e-12;
constexpr double GOLDSTONE_THRESHOLD = 1e-10;

inline int clamp_root_tol_bits(int root_tol_bits) {
    if (root_tol_bits < 1) {
        return 1;
    }
    if (root_tol_bits > MAX_EFFECTIVE_ROOT_TOL_BITS) {
        return MAX_EFFECTIVE_ROOT_TOL_BITS;
    }
    return root_tol_bits;
}

inline void require_1d_array(const py::buffer_info& buffer, const std::string& name) {
    if (buffer.ndim != 1) {
        throw py::value_error(name + " must be a one-dimensional array.");
    }
}

inline void require_non_empty_array(const py::buffer_info& buffer, const std::string& name) {
    if (buffer.size <= 0) {
        throw py::value_error(name + " must not be empty.");
    }
}

inline void require_finite_array(const double* ptr, py::ssize_t size, const std::string& name) {
    for (py::ssize_t i = 0; i < size; ++i) {
        if (!std::isfinite(ptr[i])) {
            throw py::value_error(name + " contains non-finite values.");
        }
    }
}

}  // namespace detail


struct ContinuedFractionData {
    const double a_infinity;
    const double b_infinity_squared;

    const py_array_double continuum_boundaries_squared;
    const py_array_double A;
    const py_array_double B;

    const int termination_index;
    const bool with_terminator;

    const double* A_ptr = nullptr;
    const double* B_ptr = nullptr;
    const double* continuum_boundaries_ptr = nullptr;

    py::ssize_t A_size = 0;
    py::ssize_t B_size = 0;
    py::ssize_t continuum_boundaries_size = 0;

    ContinuedFractionData(
        double a_inf,
        double b_inf_squared,
        const py_array_double& continuum_boundaries_squared_,
        const py_array_double& a,
        const py_array_double& b,
        int term_index,
        bool with_terminator_
    )
        : a_infinity(a_inf),
          b_infinity_squared(b_inf_squared),
          continuum_boundaries_squared(continuum_boundaries_squared_),
          A(a),
          B(b),
          termination_index(term_index),
          with_terminator(with_terminator_) {
        auto continuum_boundariesbuf = continuum_boundaries_squared.request();
        auto Abuf = A.request();
        auto Bbuf = B.request();

        detail::require_1d_array(continuum_boundariesbuf, "continuum_boundaries_squared");
        detail::require_1d_array(Abuf, "A");
        detail::require_1d_array(Bbuf, "B");

        detail::require_non_empty_array(Abuf, "A");
        detail::require_non_empty_array(Bbuf, "B");

        if (continuum_boundariesbuf.size != 2) {
            throw py::value_error(
                "continuum_boundaries_squared must contain exactly two values."
            );
        }

        continuum_boundaries_size = continuum_boundariesbuf.size;
        A_size = Abuf.size;
        B_size = Bbuf.size;

        continuum_boundaries_ptr = static_cast<const double*>(continuum_boundariesbuf.ptr);
        A_ptr = static_cast<const double*>(Abuf.ptr);
        B_ptr = static_cast<const double*>(Bbuf.ptr);

        detail::require_finite_array(continuum_boundaries_ptr, continuum_boundaries_size, "continuum_boundaries_squared");
        detail::require_finite_array(A_ptr, A_size, "A");
        detail::require_finite_array(B_ptr, B_size, "B");

        if (!std::isfinite(a_infinity)) {
            throw py::value_error("a_infinity must be finite.");
        }

        if (!std::isfinite(b_infinity_squared)) {
            throw py::value_error("b_infinity_squared must be finite.");
        }

        if (b_infinity_squared == 0.0) {
            throw py::value_error("b_infinity_squared must be non-zero.");
        }

        if (continuum_boundaries_ptr[0] > continuum_boundaries_ptr[1]) {
            throw py::value_error(
                "continuum_boundaries_squared must be ordered as [lower, upper]."
            );
        }

        if (B_size < A_size) {
            throw py::value_error(
                "B must contain at least as many entries as A. "
                "The continued-fraction recursion accesses B[k + 1]."
            );
        }

        validate_termination_index(termination_index);
    }

    inline double get_root(py::ssize_t idx) const {
        if (idx < 0 || idx >= continuum_boundaries_size) {
            throw py::index_error("Continuum boundary index out of range.");
        }
        return continuum_boundaries_ptr[idx];
    }

    inline void validate_termination_index(int idx) const {
        if (idx < 0 || idx >= A_size) {
            throw py::value_error(
                "termination_index is out of range. Expected "
                "0 <= termination_index < len(A)."
            );
        }

        if (idx >= B_size) {
            throw py::value_error(
                "termination_index is incompatible with len(B). "
                "Expected termination_index < len(B)."
            );
        }
    }
};


namespace detail {

std::complex<double> terminator_impl(double x, const ContinuedFractionData& data) {
    const double x_squared = x * x;
    const double p = x_squared - data.a_infinity;
    const double radicand = p * p - 4.0 * data.b_infinity_squared;

    const std::complex<double> root =
        radicand >= 0.0
            ? std::complex<double>{std::sqrt(radicand), 0.0}
            : std::complex<double>{0.0, std::sqrt(std::abs(radicand))};

    const auto boundary_index = static_cast<py::ssize_t>(x <= 0.0);

    if (x_squared > data.get_root(boundary_index)) {
        return (p - root) / (2.0 * data.b_infinity_squared);
    }

    return (p + root) / (2.0 * data.b_infinity_squared);
}


// Only the real part of x is used for the terminator.
std::vector<std::complex<double>> terminator_vector(
    const std::complex<double>* ptr_x,
    py::ssize_t n_x,
    const ContinuedFractionData& data
) {
    std::vector<std::complex<double>> result(static_cast<std::size_t>(n_x));

    for (py::ssize_t i = 0; i < n_x; ++i) {
        result[static_cast<std::size_t>(i)] = terminator_impl(ptr_x[i].real(), data);
    }

    return result;
}


std::vector<std::complex<double>> zero_terminator_vector(py::ssize_t n_x) {
    return std::vector<std::complex<double>>(
        static_cast<std::size_t>(n_x),
        std::complex<double>{0.0, 0.0}
    );
}


// This function should only be called outside the continuum.
//
// Inside the continuum the terminator is complex and the continued fraction
// generally has a finite imaginary part. In that regime, this real-valued
// denominator would not represent the full denominator correctly.
double subgap_real_denominator(
    const ContinuedFractionData& data,
    double x_squared
) {
    const double p = x_squared - data.a_infinity;

    const double terminator_value =
        data.with_terminator
            ? (p + std::sqrt(p * p - 4.0 * data.b_infinity_squared))
                  / (2.0 * data.b_infinity_squared)
            : 0.0;

    double result =
        x_squared
        - data.A_ptr[data.termination_index]
        - data.b_infinity_squared * terminator_value;

    for (int k = data.termination_index - 1; k >= 0; --k) {
        result =
            x_squared
            - data.A_ptr[k]
            - data.B_ptr[k + 1] / result;
    }

    return result;
}


double single_imaginary_part(
    const ContinuedFractionData& data,
    const std::complex<double>& x_squared
) {
    const double p = x_squared.real() - data.a_infinity;

    const double terminator_value =
        data.with_terminator
            ? (p + std::sqrt(p * p - 4.0 * data.b_infinity_squared))
                  / (2.0 * data.b_infinity_squared)
            : 0.0;

    std::complex<double> result =
        x_squared
        - data.A_ptr[data.termination_index]
        - data.b_infinity_squared * terminator_value;

    for (int k = data.termination_index - 1; k >= 0; --k) {
        result =
            x_squared
            - data.A_ptr[k]
            - data.B_ptr[k + 1] / result;
    }

    return (1.0 / result).imag();
}


void raw_denominator(
    const ContinuedFractionData& data,
    const std::complex<double>* ptr_x,
    const std::complex<double>* ptr_terminator,
    std::complex<double>* ptr_result,
    py::ssize_t n_x,
    int termination_index
) {
    data.validate_termination_index(termination_index);

    for (py::ssize_t i = 0; i < n_x; ++i) {
        const std::complex<double> x_squared = ptr_x[i] * ptr_x[i];

        ptr_result[i] =
            x_squared
            - data.A_ptr[termination_index]
            - data.b_infinity_squared * ptr_terminator[i];

        for (int k = termination_index - 1; k >= 0; --k) {
            ptr_result[i] =
                x_squared
                - data.A_ptr[k]
                - data.B_ptr[k + 1] / ptr_result[i];
        }
    }
}


void require_complex_input_1d(const py_array_cmplx& x) {
    auto xbuf = x.request();

    require_1d_array(xbuf, "x");
    require_non_empty_array(xbuf, "x");
}


}  // namespace detail


py_array_cmplx terminator(
    const py_array_double& x,
    const ContinuedFractionData& data
) {
    auto xbuf = x.request();

    detail::require_1d_array(xbuf, "x");
    detail::require_non_empty_array(xbuf, "x");

    const auto n_x = xbuf.size;

    const double* xptr = static_cast<const double*>(xbuf.ptr);

    py_array_cmplx result(n_x);
    auto rbuf = result.request();

    auto* rptr = static_cast<std::complex<double>*>(rbuf.ptr);

    {
        py::gil_scoped_release release;

        for (py::ssize_t i = 0; i < n_x; ++i) {
            rptr[i] = detail::terminator_impl(xptr[i], data);
        }
    }

    return result;
}


py_array_cmplx denominator(
    const py_array_cmplx& x,
    const ContinuedFractionData& data
) {
    auto xbuf = x.request();

    detail::require_1d_array(xbuf, "x");
    detail::require_non_empty_array(xbuf, "x");

    if (data.B_ptr[0] == 0.0) {
        throw py::value_error("B[0] must be non-zero for denominator normalization.");
    }

    const auto n_x = xbuf.size;

    const auto* ptr_x =
        static_cast<const std::complex<double>*>(xbuf.ptr);

    py_array_cmplx result(n_x);
    auto rbuf = result.request();

    auto* ptr_result =
        static_cast<std::complex<double>*>(rbuf.ptr);

    {
        py::gil_scoped_release release;

        const auto terminator_data =
            data.with_terminator
                ? detail::terminator_vector(ptr_x, n_x, data)
                : detail::zero_terminator_vector(n_x);

        detail::raw_denominator(
            data,
            ptr_x,
            terminator_data.data(),
            ptr_result,
            n_x,
            data.termination_index
        );

        for (py::ssize_t i = 0; i < n_x; ++i) {
            ptr_result[i] /= data.B_ptr[0];
        }
    }

    return result;
}


py_array_cmplx continued_fraction(
    const py_array_cmplx& x,
    const ContinuedFractionData& data
) {
    auto xbuf = x.request();

    detail::require_1d_array(xbuf, "x");
    detail::require_non_empty_array(xbuf, "x");

    const auto n_x = xbuf.size;

    const auto* ptr_x =
        static_cast<const std::complex<double>*>(xbuf.ptr);

    py_array_cmplx result(n_x);
    auto rbuf = result.request();

    auto* ptr_result =
        static_cast<std::complex<double>*>(rbuf.ptr);

    {
        py::gil_scoped_release release;

        const auto terminator_data =
            data.with_terminator
                ? detail::terminator_vector(ptr_x, n_x, data)
                : detail::zero_terminator_vector(n_x);

        detail::raw_denominator(
            data,
            ptr_x,
            terminator_data.data(),
            ptr_result,
            n_x,
            data.termination_index
        );

        for (py::ssize_t i = 0; i < n_x; ++i) {
            ptr_result[i] = data.B_ptr[0] / ptr_result[i];
        }
    }

    return result;
}


py_array_cmplx continued_fraction_varied_depth(
    const py_array_cmplx& x,
    const ContinuedFractionData& data,
    const py::array_t<int, pyarray_flags>& shift_range
) {
    auto xbuf = x.request();
    auto shiftbuf = shift_range.request();

    detail::require_1d_array(xbuf, "x");
    detail::require_non_empty_array(xbuf, "x");

    detail::require_1d_array(shiftbuf, "shift_range");
    detail::require_non_empty_array(shiftbuf, "shift_range");

    const auto n_x = xbuf.size;
    const auto n_shift_range = shiftbuf.size;

    const auto* ptr_x =
        static_cast<const std::complex<double>*>(xbuf.ptr);

    const auto* ptr_shift_range =
        static_cast<const int*>(shiftbuf.ptr);

    for (py::ssize_t r = 0; r < n_shift_range; ++r) {
        const int current_termination_index =
            data.termination_index + ptr_shift_range[r];

        data.validate_termination_index(current_termination_index);
    }

    py_array_cmplx result({n_shift_range, n_x});
    auto rbuf = result.request();

    auto* ptr_result =
        static_cast<std::complex<double>*>(rbuf.ptr);

    {
        py::gil_scoped_release release;

        const auto terminator_data =
            data.with_terminator
                ? detail::terminator_vector(ptr_x, n_x, data)
                : detail::zero_terminator_vector(n_x);

        std::vector<std::complex<double>> x_squared_data(
            static_cast<std::size_t>(n_x)
        );

        for (py::ssize_t i = 0; i < n_x; ++i) {
            x_squared_data[static_cast<std::size_t>(i)] = ptr_x[i] * ptr_x[i];
        }

        for (py::ssize_t r = 0; r < n_shift_range; ++r) {
            const int current_termination_index =
                data.termination_index + ptr_shift_range[r];

            for (py::ssize_t i = 0; i < n_x; ++i) {
                const auto x_squared =
                    x_squared_data[static_cast<std::size_t>(i)];

                auto& current_result = ptr_result[r * n_x + i];

                current_result =
                    x_squared
                    - data.A_ptr[current_termination_index]
                    - data.b_infinity_squared * terminator_data[static_cast<std::size_t>(i)];

                for (int k = current_termination_index - 1; k >= 0; --k) {
                    current_result =
                        x_squared
                        - data.A_ptr[k]
                        - data.B_ptr[k + 1] / current_result;
                }

                current_result = data.B_ptr[0] / current_result;
            }
        }
    }

    return result;
}


// This classification is designed for finite-energy simple poles.
// It is known to be less reliable for Goldstone peaks / double-pole-like
// structures near omega = 0.
std::vector<std::pair<double, double>> classify_bound_states(
    const ContinuedFractionData& data,
    std::size_t n_scan,
    double weight_domega,
    int root_tol_bits,
    std::uintmax_t max_iter
) {
    if (n_scan < 3) {
        throw py::value_error("n_scan must be at least 3.");
    }

    if (weight_domega <= 0.0 || !std::isfinite(weight_domega)) {
        throw py::value_error("weight_domega must be positive and finite.");
    }

    root_tol_bits = detail::clamp_root_tol_bits(root_tol_bits);

    const double root_tol = std::max(
        std::ldexp(1.0, 1 - root_tol_bits),
        4.0 * std::numeric_limits<double>::epsilon()
    );

    const double dz =
        detail::BOOST_SHIFT
        * (data.get_root(0) - weight_domega - root_tol)
        / static_cast<double>(n_scan);

    std::vector<std::pair<double, double>> results;

    if (dz <= 0.0) {
        return results;
    }

    auto denom = [&data](double z) -> double {
        return detail::subgap_real_denominator(data, z);
    };

    std::complex<double> z_squared = {0.0, detail::EPS};

    auto minimize_function =
        [&data, &z_squared](double z) -> double {
            z_squared.real(z);
            const double f = detail::single_imaginary_part(data, z_squared);
            return std::atan(f);
        };

    auto set_coefficient_of_pole =
        [&]() {
            auto& last_result = results.back();

            const double peak_position = last_result.first;

            // For finite-energy simple poles, the spectral weight is the
            // residue. It is approximated by a symmetric finite difference.
            if (peak_position > detail::GOLDSTONE_THRESHOLD + root_tol) {
                const double omega_plus = peak_position + weight_domega;
                const double omega_minus = peak_position - weight_domega;

                last_result.second =
                    2.0
                    * weight_domega
                    * data.B_ptr[0]
                    / (
                        denom(omega_plus * omega_plus)
                        - denom(omega_minus * omega_minus)
                    );
            } else {
                // For an omega = 0 mode with G(omega) ~ Z / omega^2,
                // use a finite-difference form inspired by l'Hopital's rule:
                //
                //     Z = lim_{omega -> 0} omega^2 / denom(omega^2).
                //
                // Numerically this remains delicate for Goldstone-like peaks.
                last_result.second =
                    weight_domega
                    * weight_domega
                    * data.B_ptr[0]
                    / (
                        4.0
                        * denom(weight_domega * weight_domega)
                    );
            }
        };

    double f_left = minimize_function(0.0);
    double f_center = minimize_function(dz);
    double f_right = 0.0;

    {
        py::gil_scoped_release release;

        for (std::size_t i = 2; i < n_scan; ++i) {
            f_right = minimize_function(static_cast<double>(i) * dz);

            if (f_left > f_center && f_right > f_center) {
                // There is a local minimum somewhere between the neighboring
                // scan points.
                const auto minimization_result =
                    boost::math::tools::brent_find_minima(
                        minimize_function,
                        static_cast<double>(i - 2) * dz,
                        static_cast<double>(i) * dz,
                        root_tol_bits
                    );

                double search_range = root_tol;

                double a = minimization_result.first - search_range;
                double b = minimization_result.first + search_range;

                double fa = denom(a);
                double fb = denom(b);

                while (
                    std::signbit(fa) == std::signbit(fb)
                    && search_range > detail::MIN_SEARCH_RANGE
                ) {
                    search_range /= 10.0;
                    a = minimization_result.first - search_range;
                    fa = denom(a);
                }

                search_range = root_tol;

                while (
                    std::signbit(fa) == std::signbit(fb)
                    && search_range > detail::MIN_SEARCH_RANGE
                ) {
                    search_range /= 10.0;
                    b = minimization_result.first + search_range;
                    fb = denom(b);
                }

                if (search_range > detail::MIN_SEARCH_RANGE) {
                    std::uintmax_t iter = max_iter;

                    const auto root_result =
                        boost::math::tools::toms748_solve(
                            denom,
                            a,
                            b,
                            boost::math::tools::eps_tolerance<double>(
                                root_tol_bits
                            ),
                            iter
                        );

                    const double root_z =
                        0.5 * (root_result.first + root_result.second);

                    if (root_z >= 0.0) {
                        results.emplace_back(std::sqrt(root_z), 0.0);
                        set_coefficient_of_pole();
                    }
                }
            }

            f_left = f_center;
            f_center = f_right;
        }
    }

    return results;
}


PYBIND11_MODULE(cpp_continued_fraction, m) {
    m.doc() = "Fast continued-fraction routines implemented in C++ with pybind11.";

    py::class_<ContinuedFractionData>(m, "ContinuedFractionData")
        .def(
            py::init<
                double,
                double,
                const py_array_double&,
                const py_array_double&,
                const py_array_double&,
                int,
                bool
            >(),
            py::arg("a_infinity"),
            py::arg("b_infinity_squared"),
            py::arg("continuum_boundaries_squared"),
            py::arg("A"),
            py::arg("B"),
            py::arg("termination_index"),
            py::arg("with_terminator"),
            R"pbdoc(
                Data container for continued-fraction evaluation.

                Parameters
                ----------
                a_infinity : float
                    Asymptotic value of the A coefficients.
                b_infinity_squared : float
                    Square of the asymptotic B coefficient.
                continuum_boundaries_squared : array_like, shape (2,)
                    Squared lower and upper continuum boundaries.
                A : array_like
                    Continued-fraction A coefficients.
                B : array_like
                    Continued-fraction B coefficients.
                termination_index : int
                    Index at which the continued fraction is terminated.
                with_terminator : bool
                    Whether to attach the square-root terminator.
            )pbdoc"
        )
        .def_readonly("a_infinity", &ContinuedFractionData::a_infinity)
        .def_readonly("b_infinity_squared", &ContinuedFractionData::b_infinity_squared)
        .def_readonly(
            "continuum_boundaries_squared",
            &ContinuedFractionData::continuum_boundaries_squared
        )
        .def_readonly("A", &ContinuedFractionData::A)
        .def_readonly("B", &ContinuedFractionData::B)
        .def_readonly("termination_index", &ContinuedFractionData::termination_index)
        .def_readonly("with_terminator", &ContinuedFractionData::with_terminator)
        .def_readonly("A_size", &ContinuedFractionData::A_size)
        .def_readonly("B_size", &ContinuedFractionData::B_size);

    m.def(
        "terminator",
        &terminator,
        py::arg("x"),
        py::arg("data"),
        R"pbdoc(
            Compute the square-root terminator.

            Parameters
            ----------
            x : array_like
                Real frequencies.
            data : ContinuedFractionData
                Continued-fraction data.

            Returns
            -------
            numpy.ndarray
                Complex terminator values.
        )pbdoc"
    );

    m.def(
        "denominator",
        &denominator,
        py::arg("x"),
        py::arg("data"),
        R"pbdoc(
            Compute the normalized denominator of the continued fraction.

            Parameters
            ----------
            x : array_like
                Complex frequencies.
            data : ContinuedFractionData
                Continued-fraction data.

            Returns
            -------
            numpy.ndarray
                Complex denominator values divided by B[0].
        )pbdoc"
    );

    m.def(
        "continued_fraction",
        &continued_fraction,
        py::arg("x"),
        py::arg("data"),
        R"pbdoc(
            Compute the continued fraction.

            Parameters
            ----------
            x : array_like
                Complex frequencies.
            data : ContinuedFractionData
                Continued-fraction data.

            Returns
            -------
            numpy.ndarray
                Complex continued-fraction values.
        )pbdoc"
    );

    m.def(
        "continued_fraction_varied_depth",
        &continued_fraction_varied_depth,
        py::arg("x"),
        py::arg("data"),
        py::arg("shift_range"),
        R"pbdoc(
            Compute the continued fraction for several shifted termination depths.

            Parameters
            ----------
            x : array_like
                Complex frequencies.
            data : ContinuedFractionData
                Continued-fraction data.
            shift_range : array_like of int
                Integer shifts added to data.termination_index.

            Returns
            -------
            numpy.ndarray
                Complex array of shape ``(len(shift_range), len(x))``.
        )pbdoc"
    );

    m.def(
        "classify_bound_states",
        &classify_bound_states,
        py::arg("data"),
        py::arg("n_scan"),
        py::arg("weight_domega"),
        py::arg("root_tol_bits"),
        py::arg("max_iter"),
        R"pbdoc(
            Compute approximate energies and spectral weights of bound states.

            This routine scans the interval below the lower continuum edge and
            searches for zeros of the real denominator.

            Notes
            -----
            The method assumes that finite-energy bound states are simple poles.
            Goldstone-like peaks near omega = 0 are numerically delicate.

            Returns
            -------
            list[tuple[float, float]]
                Pairs ``(energy, weight)``.
        )pbdoc"
    );
}
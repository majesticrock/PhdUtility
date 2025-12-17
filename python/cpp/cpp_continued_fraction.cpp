#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <complex>
#include <vector>
#include <cmath>
#include <list>
#include <limits>
#include <boost/math/tools/roots.hpp>

#define MROCK_CF_INIT(result_dims) py::buffer_info buf_x = x.request(); \
    std::complex<double>* __restrict ptr_x = static_cast<std::complex<double>*>(buf_x.ptr); \
    const size_t n_x = buf_x.size; \
    std::vector<std::complex<double>> terminator_data = data.with_terminator ? detail::terminator_vector(ptr_x, n_x, data) : std::vector<std::complex<double>>(n_x, std::complex<double>{0.,0.}); \
    \
    py_array_cmplx result( result_dims ); \
    py::buffer_info buf_result = result.request(); \
    std::complex<double>* __restrict ptr_result = static_cast<std::complex<double>*>(buf_result.ptr); 

#define M_KOMMA ,

namespace py = pybind11;
constexpr auto pyarray_flags = py::array::c_style | py::array::forcecast;
using py_array_double = py::array_t<double, pyarray_flags>;
using py_array_cmplx  = py::array_t<std::complex<double>, pyarray_flags>;

struct ContinuedFractionData {
    const double a_infinity;
    const double b_infinity_squared;
    const py_array_double continuum_boundaries_squared;
    const py_array_double A;
    const py_array_double B;

    const double* __restrict A_ptr;
    const double* __restrict B_ptr;
    const double* __restrict cb_ptr;

    const int termination_index;
    const bool with_terminator;

    ContinuedFractionData(double a_inf, double b_inf_squared, 
        const py_array_double& r, const py_array_double& a, const py_array_double& b, 
        int term_index, bool _with_terminator)
        : a_infinity(a_inf), 
        b_infinity_squared(b_inf_squared), 
        continuum_boundaries_squared(r),
        A(a), 
        B(b), 
        termination_index(term_index),
        with_terminator(_with_terminator)
    {
        auto cbbuf = continuum_boundaries_squared.request();
        auto Abuf = A.request();
        auto Bbuf = B.request();
        cb_ptr = static_cast<const double*>(cbbuf.ptr);
        A_ptr    = static_cast<const double*>(Abuf.ptr);
        B_ptr    = static_cast<const double*>(Bbuf.ptr);
    }

    inline double getRoot(int idx) const {
        return cb_ptr[idx];
    }
};

namespace detail {
    inline std::complex<double> terminator(const double x, ContinuedFractionData const& data) {
        const double x_squared{ x * x };
        const double p{ x_squared - data.a_infinity };

        const double radicant{p * p - 4. * data.b_infinity_squared};
        const std::complex<double> root = radicant >= 0.0 ? std::complex<double>{sqrt(radicant), 0.0} : std::complex<double>{0.0, sqrt(abs(radicant))};

        if (x_squared > data.getRoot(static_cast<int>(x <= 0))) {
            return (p - root) / (2. * data.b_infinity_squared);
        }
        else {
            return (p + root) / (2. * data.b_infinity_squared);
        }
    }

    // Only the real part of x is considered
    inline std::vector<std::complex<double>> terminator_vector(const std::complex<double>* ptr_x, const size_t n_x, ContinuedFractionData const& data)
    {
        std::vector<std::complex<double>> ret(n_x);
        for (size_t i = 0U; i < n_x; ++i) {
            const double x_squared{ ptr_x[i].real() * ptr_x[i].real() };
            const double p{ x_squared - data.a_infinity };

            const double radicant{p * p - 4. * data.b_infinity_squared};
            const std::complex<double> root = radicant >= 0.0 ? std::complex<double>{sqrt(radicant), 0.0} : std::complex<double>{0.0, sqrt(abs(radicant))};

            if (x_squared > data.getRoot(static_cast<int>(ptr_x[i].real() <= 0))) {
                ret[i] = (p - root) / (2. * data.b_infinity_squared);
            }
            else {
                ret[i] = (p + root) / (2. * data.b_infinity_squared);
            }
        }
        return ret;
    }

    // This function may only be called outside of the continuum
    // Inside the continuum the terminator will be a complex number
    // and consequently, the continued fraction will have a finite imaginary part,
    // which would make this function return wrong values
    inline double subgap_real_denominator(ContinuedFractionData const& data, const double x_squared)
    {
        const double p{ x_squared - data.a_infinity };
        const double terminator_value = (p + sqrt(p * p - 4. * data.b_infinity_squared)) / (2. * data.b_infinity_squared);

        double ret = x_squared - data.A_ptr[data.termination_index] - data.B_ptr[data.termination_index + 1] * terminator_value;
        for (int k = data.termination_index - 1; k >= 0; --k) {
            ret = x_squared - data.A_ptr[k] - data.B_ptr[k+1] / ret;
        }
        return ret;
    }

    inline void raw_denominator(ContinuedFractionData const& data, 
        const std::complex<double>* __restrict ptr_x, const std::complex<double>* __restrict ptr_terminator, 
        std::complex<double>* __restrict ptr_result, const size_t n_x)
    {
        for (size_t i = 0U; i < n_x; ++i) {
            const std::complex<double> x_squared{ ptr_x[i] * ptr_x[i] };
            ptr_result[i] = x_squared - data.A_ptr[data.termination_index] - data.B_ptr[data.termination_index + 1] * ptr_terminator[i];
            for (int k = data.termination_index - 1; k >= 0; --k) {
                ptr_result[i] = x_squared - data.A_ptr[k] - data.B_ptr[k+1] / ptr_result[i];
            }
        }
    }
}

py_array_cmplx denominator(const py_array_cmplx& x, ContinuedFractionData const& data)
{
    MROCK_CF_INIT(n_x);
    detail::raw_denominator(data, ptr_x, terminator_data.data(), ptr_result, n_x);

    for (size_t i = 0U; i < n_x; ++i) {
        ptr_result[i] /= data.B_ptr[0];
    }
    return result;
}

py_array_cmplx continued_fraction(const py_array_cmplx& x, ContinuedFractionData const& data)
{
    MROCK_CF_INIT(n_x);
    detail::raw_denominator(data, ptr_x, terminator_data.data(), ptr_result, n_x);

    for (size_t i = 0U; i < n_x; ++i) {
        ptr_result[i] = data.B_ptr[0] / ptr_result[i];
    }
    return result;
}

py_array_cmplx continued_fraction_varied_depth(const py_array_cmplx& x, 
    ContinuedFractionData const& data, const py::array_t<int, pyarray_flags> shift_range)
{
    py::buffer_info buf_shift_range = shift_range.request();
    int* __restrict ptr_shift_range = static_cast<int*>(buf_shift_range.ptr);
    const size_t n_shift_range = buf_shift_range.size;

    MROCK_CF_INIT({n_shift_range M_KOMMA n_x});

    std::vector<std::complex<double>> x_squared_data(n_x);
    for (size_t i = 0U; i < n_x; ++i) {
        x_squared_data[i]  = ptr_x[i] * ptr_x[i];
    }

    for (size_t r = 0U; r < n_shift_range; ++r) {
        const int current_termination_index = data.termination_index + ptr_shift_range[r];
        for (size_t i = 0U; i < n_x; ++i) {
            ptr_result[r * n_x + i] = (ptr_x[i] * ptr_x[i]) - data.A_ptr[current_termination_index] - data.B_ptr[current_termination_index + 1] * terminator_data[i];
            for (int k = current_termination_index - 1; k >= 0; --k) {
                ptr_result[r * n_x + i] = x_squared_data[i] - data.A_ptr[k] - data.B_ptr[k+1] / ptr_result[r * n_x + i];
            }
            ptr_result[r * n_x + i] = data.B_ptr[0] / ptr_result[r * n_x + i];
        }
    }

    return result;
}

std::list<std::pair<double, double>> classify_bound_states(ContinuedFractionData const& data, size_t n_scan, double weight_domega, int root_tol_bits, int maxiter)
{
    const double root_tol = std::max(ldexp(1.0, 1-root_tol_bits), 4 * std::numeric_limits<double>::epsilon());
    const double dz = (data.getRoot(0) - root_tol) / n_scan;
    double z_sqr{};

    auto denom = [&data](double z) -> double {
        return detail::subgap_real_denominator(data, z);
    };

    double a = denom(z_sqr);
    double b;

    // saves the pairs {peak position, weight}
    std::list<std::pair<double, double>> results;

    for (size_t i = 0U; i < n_scan; ++i) {
        z_sqr = (i + 1U) * dz;
        b = denom(z_sqr);

        if (std::signbit(a) != std::signbit(b)) {
            // Root in interval [z_sqr - dz, z_sqr]
            results.emplace_back(std::pair<double, double>{0., 0.});

            std::uintmax_t boost_max_it{100U};
            const auto sol = boost::math::tools::toms748_solve(denom, z_sqr - dz, z_sqr, boost::math::tools::eps_tolerance<double>(root_tol_bits), boost_max_it);
            results.back().first = sqrt(0.5 * (sol.first + sol.second));

            // The spectral weight of the bound state is its residue
            // Besides the Goldstone bosons, all peaks (thus far) are delta peaks
            // i.e., the real part is a simple pole.
            // Thus, the residue is given by the below expression
            if (results.back().first > (1e-10 + root_tol)) {
                results.back().second = 2 * weight_domega * data.B_ptr[0] / (
                          denom((results.back().first + weight_domega) * (results.back().first + weight_domega)) 
                        - denom((results.back().first - weight_domega) * (results.back().first - weight_domega)) 
                    );
            }
            else {
                // The residue of 1/omega^2 is 0, thus we revert to l'Hopital:
                // lim_w->0 w^2 / denom(w^2) = Z
                //     = 2 / (d^2 denom(w^2) / dw^2)
                results.back().second = weight_domega * weight_domega * data.B_ptr[0] / (2. * 
                    2. * denom(weight_domega * weight_domega) // denom((z_0 + dz)^2) + denom((z_0 - dz)^2) with z_0 = 0
                    // - 2 * denom(0) = 0
                );
            }
        }
        a = b;
    }

    return results;
}


PYBIND11_MODULE(cpp_continued_fraction, m) {
    // Binding for ContinuedFractionData class
    py::class_<ContinuedFractionData>(m, "ContinuedFractionData")
        .def(py::init<double, double, const py_array_double&, const py_array_double&, const py_array_double&, int, bool>())
        .def_readonly("a_infinity", &ContinuedFractionData::a_infinity)
        .def_readonly("b_infinity_squared", &ContinuedFractionData::b_infinity_squared)
        .def_readonly("continuum_boundaries_squared", &ContinuedFractionData::continuum_boundaries_squared)
        .def_readonly("A", &ContinuedFractionData::A)
        .def_readonly("B", &ContinuedFractionData::B)
        .def_readonly("termination_index", &ContinuedFractionData::termination_index)
        .def_readonly("with_terminator", &ContinuedFractionData::with_terminator);

    // Binding of functions
    m.def("denominator", &continued_fraction, 
        "Compute the denominator of the continued fraction based on input complex numbers and provided data.",
        py::arg("x"),
        py::arg("data")
    );

    m.def("continued_fraction", &continued_fraction, 
        "Compute the continued fraction based on input complex numbers and provided data.",
        py::arg("x"),
        py::arg("data")
    );

    m.def("continued_fraction_varied_depth", &continued_fraction_varied_depth, 
        "Compute the continued fraction based on input complex numbers and provided data."
        "The termination depth is varied by the indizes provided in shift_range."
        "Returns a 2D array of shape (R, N), where N is len(x) and R len(shift_range).",
        py::arg("x"),
        py::arg("data"),
        py::arg("shift_range")
    );

    m.def("classify_bound_states", &classify_bound_states,
        "Computes the energies and spectral weights of bound states, i.e., states with energy [0, omega_-]."
        "Assumes that any bound state with finite energy is a delta distribution -> simple pole"
        "and that any bound state with omega_0 = 0 (e.g., the Goldstones) is the derivative of a delta distribution"
        "-> G(omega) ~ 1/omega^2."
        "\n"
        "Returns a list of pairs. The first element of the pairs is the energy."
        "The second element is the weight (or in case of delta' the prefactor).",
        py::arg("data"),
        py::arg("n_scan"),
        py::arg("weight_domega"),
        py::arg("root_tol_bits"),
        py::arg("maxiter")
    );
}

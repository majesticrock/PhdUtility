#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <complex>
#include <vector>
#include <cmath>

namespace py = pybind11;

struct ContinuedFractionData {
    const double a_infinity;
    const double b_infinity_squared;
    const py::array_t<double> roots;
    const py::array_t<double> A;
    const py::array_t<double> B;
    const int termination_index;

    ContinuedFractionData(double a_inf, double b_inf_squared, const py::array_t<double>& r, 
        const py::array_t<double>& a, const py::array_t<double>& b, int term_index)
        : a_infinity(a_inf), b_infinity_squared(b_inf_squared), roots(r), A(a), B(b), termination_index(term_index) {}

    inline double getRoot(int idx) const {
        return roots.unchecked<1>()(idx);
    }
};

std::complex<double> terminator(const double x, ContinuedFractionData const& data) {
    auto root = [&data](double _p) -> std::complex<double> {
        return sqrt(std::complex<double>(_p * _p - 4. * data.b_infinity_squared));
    };

    const double x_squared{ x * x };
    const double p{ x_squared - data.a_infinity };

    if(x_squared > data.getRoot(static_cast<int>(x <= 0))) {
        return (p - root(p)) / (2. * data.b_infinity_squared);
    }
    else {
        return (p + root(p)) / (2. * data.b_infinity_squared);
    }
}

py::array_t<std::complex<double>> continued_fraction(const py::array_t<std::complex<double>> x, ContinuedFractionData const& data, const bool with_terminator)
{
    py::buffer_info buf_x = x.request();
    std::complex<double>* ptr_x = static_cast<std::complex<double>*>(buf_x.ptr);
    const size_t n_x = buf_x.size;

    auto result = py::array_t<std::complex<double>>(n_x);
    auto ptr_result = result.mutable_unchecked<1>();

    const auto ptr_A = data.A.unchecked<1>();
    const auto ptr_B = data.B.unchecked<1>();

    for(size_t i = 0U; i < n_x; ++i) {
        const std::complex<double> x_squared{ ptr_x[i] * ptr_x[i] };
        ptr_result(i) = (ptr_x[i] * ptr_x[i]) - ptr_A(data.termination_index)
                        - (with_terminator ? ptr_B(data.termination_index + 1) * terminator(ptr_x[i].real(), data) : double{});
        for(int k = data.termination_index - 1; k >= 0; --k) {
            ptr_result(i) = x_squared - ptr_A(k) - ptr_B(k+1) / ptr_result(i);
        }
        ptr_result(i) = ptr_B(0) / ptr_result(i);
    }

    return result;
}

PYBIND11_MODULE(cpp_continued_fraction, m) {
    // Binding for ContinuedFractionData class
    py::class_<ContinuedFractionData>(m, "ContinuedFractionData")
        .def(py::init<double, double, const py::array_t<double>&, const py::array_t<double>&, const py::array_t<double>&, int>())
        .def_readonly("a_infinity", &ContinuedFractionData::a_infinity)
        .def_readonly("b_infinity_squared", &ContinuedFractionData::b_infinity_squared)
        .def_readonly("roots", &ContinuedFractionData::roots)
        .def_readonly("A", &ContinuedFractionData::A)
        .def_readonly("B", &ContinuedFractionData::B)
        .def_readonly("termination_index", &ContinuedFractionData::termination_index);

    // Binding for the continued_fraction function
    m.def("continued_fraction", &continued_fraction, 
        "Compute the continued fraction based on input complex numbers and provided data.",
        py::arg("x"),
        py::arg("data"),
        py::arg("with_terminator") = false);
}

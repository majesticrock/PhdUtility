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

    const double* A_ptr = nullptr;
    const double* B_ptr = nullptr;
    const double* root_ptr = nullptr;

    ContinuedFractionData(double a_inf, double b_inf_squared, const py::array_t<double>& r, const py::array_t<double>& a, const py::array_t<double>& b, int term_index)
        : a_infinity(a_inf), 
        b_infinity_squared(b_inf_squared), 
        roots(r),
        A(a), 
        B(b), 
        termination_index(term_index) 
    {
        auto rbuf = roots.request();
        auto Abuf = A.request();
        auto Bbuf = B.request();
        root_ptr = static_cast<const double*>(rbuf.ptr);
        A_ptr    = static_cast<const double*>(Abuf.ptr);
        B_ptr    = static_cast<const double*>(Bbuf.ptr);
    }

    inline double getRoot(int idx) const {
        return root_ptr[idx];
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

    py::array_t<std::complex<double>> result(n_x);
    auto ptr_result = result.mutable_unchecked<1>();

    for (size_t i = 0U; i < n_x; ++i) {
        const std::complex<double> x_squared{ ptr_x[i] * ptr_x[i] };
        ptr_result(i) = (ptr_x[i] * ptr_x[i]) - data.A_ptr[data.termination_index]
                        - (with_terminator ? data.B_ptr[data.termination_index + 1] * terminator(ptr_x[i].real(), data) : double{});
        for (int k = data.termination_index - 1; k >= 0; --k) {
            ptr_result(i) = x_squared - data.A_ptr[k] - data.B_ptr[k+1] / ptr_result(i);
        }
        ptr_result(i) = data.B_ptr[0] / ptr_result(i);
    }

    return result;
}

py::array_t<std::complex<double>> continued_fraction_varied_depth(const py::array_t<std::complex<double>> x, 
    ContinuedFractionData const& data, const py::array_t<int> shift_range, const bool with_terminator)
{
    py::buffer_info buf_x = x.request();
    std::complex<double>* ptr_x = static_cast<std::complex<double>*>(buf_x.ptr);
    const size_t n_x = buf_x.size;

    py::buffer_info buf_shift_range = shift_range.request();
    int* ptr_shift_range = static_cast<int*>(buf_shift_range.ptr);
    const size_t n_shift_range = buf_shift_range.size;

    py::array_t<std::complex<double>> result({ n_shift_range, n_x });
    auto ptr_result = result.mutable_unchecked<2>();

    std::vector<std::complex<double>> terminator_data(n_x, std::complex<double>{0., 0.});
    if (with_terminator) {
        for (size_t i = 0U; i < n_x; ++i) 
            terminator_data[i] = terminator(ptr_x[i].real(), data);
    }
    std::vector<std::complex<double>> x_squared_data(n_x);
    for (size_t i = 0U; i < n_x; ++i) {
        x_squared_data[i]  = ptr_x[i] * ptr_x[i];
    }

    for (size_t r = 0U; r < n_shift_range; ++r) {
        const int current_termination_index = data.termination_index + ptr_shift_range[r];
        for (size_t i = 0U; i < n_x; ++i) {
            ptr_result(r, i) = (ptr_x[i] * ptr_x[i]) - data.A_ptr[current_termination_index] - data.B_ptr[current_termination_index + 1] * terminator_data[i];
            for (int k = current_termination_index - 1; k >= 0; --k) {
                ptr_result(r, i) = x_squared_data[i] - data.A_ptr[k] - data.B_ptr[k+1] / ptr_result(r, i);
            }
            ptr_result(r, i) = data.B_ptr[0] / ptr_result(r, i);
        }
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
        py::arg("with_terminator") = true);

    m.def("continued_fraction_varied_depth", &continued_fraction_varied_depth, 
        "Compute the continued fraction based on input complex numbers and provided data."
        "The termination depth is varied by the indizes provided in shift_range."
        "Returns a 2D array of shape (R, N), where N is len(x) and R len(shift_range).",
        py::arg("x"),
        py::arg("data"),
        py::arg("shift_range"),
        py::arg("with_terminator") = true);
}

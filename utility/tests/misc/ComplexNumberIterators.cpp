#include <algorithm>
#include <complex>
#include <compare>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>

#include <mrock/utility/ComplexNumberIterators.hpp>

#define FAIL_IF_NOT(cond)                                      \
    if (!(cond)) {                                             \
        std::cerr << "FAILED: " #cond                          \
                  << " at line " << __LINE__ << std::endl;     \
        return 1;                                              \
    }

int main() {
    using namespace mrock::utility;

    using Complex = std::complex<double>;

    std::vector<Complex> data{
        {1.0, 10.0},
        {2.0, 20.0},
        {3.0, 30.0},
        {4.0, 40.0}
    };

    {
        std::cout << "Test 1: real-part iterator from vector\n";

        auto begin = make_real_part_iterator(data);
        auto end = make_real_part_iterator_end(data);

        FAIL_IF_NOT(end - begin == static_cast<std::ptrdiff_t>(data.size()));

        FAIL_IF_NOT(*begin == 1.0);
        FAIL_IF_NOT(begin[0] == 1.0);
        FAIL_IF_NOT(begin[1] == 2.0);
        FAIL_IF_NOT(begin[2] == 3.0);
        FAIL_IF_NOT(begin[3] == 4.0);

        std::vector<double> real_parts(begin, end);

        FAIL_IF_NOT(real_parts.size() == data.size());
        FAIL_IF_NOT(real_parts[0] == 1.0);
        FAIL_IF_NOT(real_parts[1] == 2.0);
        FAIL_IF_NOT(real_parts[2] == 3.0);
        FAIL_IF_NOT(real_parts[3] == 4.0);

        std::cout << "Real parts:";
        for (double x : real_parts) {
            std::cout << ' ' << x;
        }
        std::cout << "\n\n";
    }

    {
        std::cout << "Test 2: imag-part iterator from vector\n";

        auto begin = make_imag_part_iterator(data);
        auto end = make_imag_part_iterator_end(data);

        FAIL_IF_NOT(end - begin == static_cast<std::ptrdiff_t>(data.size()));

        FAIL_IF_NOT(*begin == 10.0);
        FAIL_IF_NOT(begin[0] == 10.0);
        FAIL_IF_NOT(begin[1] == 20.0);
        FAIL_IF_NOT(begin[2] == 30.0);
        FAIL_IF_NOT(begin[3] == 40.0);

        std::vector<double> imag_parts(begin, end);

        FAIL_IF_NOT(imag_parts.size() == data.size());
        FAIL_IF_NOT(imag_parts[0] == 10.0);
        FAIL_IF_NOT(imag_parts[1] == 20.0);
        FAIL_IF_NOT(imag_parts[2] == 30.0);
        FAIL_IF_NOT(imag_parts[3] == 40.0);

        std::cout << "Imag parts:";
        for (double x : imag_parts) {
            std::cout << ' ' << x;
        }
        std::cout << "\n\n";
    }

    {
        std::cout << "Test 3: iterator increment and decrement\n";

        auto it = make_real_part_iterator(data);

        FAIL_IF_NOT(*it == 1.0);

        ++it;
        FAIL_IF_NOT(*it == 2.0);

        it++;
        FAIL_IF_NOT(*it == 3.0);

        --it;
        FAIL_IF_NOT(*it == 2.0);

        it--;
        FAIL_IF_NOT(*it == 1.0);

        std::cout << "Increment/decrement test passed.\n\n";
    }

    {
        std::cout << "Test 4: random-access operations\n";

        auto begin = make_real_part_iterator(data);
        auto end = make_real_part_iterator_end(data);

        auto it = begin + 2;
        FAIL_IF_NOT(*it == 3.0);

        it -= 1;
        FAIL_IF_NOT(*it == 2.0);

        it += 2;
        FAIL_IF_NOT(*it == 4.0);

        auto previous = it - 1;
        FAIL_IF_NOT(*previous == 3.0);

        FAIL_IF_NOT(end - begin == 4);
        FAIL_IF_NOT(it - begin == 3);

        FAIL_IF_NOT(begin < end);
        FAIL_IF_NOT(begin <= end);
        FAIL_IF_NOT(end > begin);
        FAIL_IF_NOT(end >= begin);
        FAIL_IF_NOT(begin == begin);
        FAIL_IF_NOT(begin != end);

        std::cout << "Random-access operation test passed.\n\n";
    }

    {
        std::cout << "Test 5: iterator construction from raw pointer\n";

        auto real_begin = make_real_part_iterator(data.data());
        auto real_end = make_real_part_iterator_end(data.data(), data.size());

        auto imag_begin = make_imag_part_iterator(data.data());
        auto imag_end = make_imag_part_iterator_end(data.data(), data.size());

        FAIL_IF_NOT(real_end - real_begin == static_cast<std::ptrdiff_t>(data.size()));
        FAIL_IF_NOT(imag_end - imag_begin == static_cast<std::ptrdiff_t>(data.size()));

        FAIL_IF_NOT(real_begin[0] == 1.0);
        FAIL_IF_NOT(real_begin[3] == 4.0);

        FAIL_IF_NOT(imag_begin[0] == 10.0);
        FAIL_IF_NOT(imag_begin[3] == 40.0);

        std::cout << "Raw pointer construction test passed.\n\n";
    }

    {
        std::cout << "Test 6: use with standard algorithms\n";

        auto real_begin = make_real_part_iterator(data);
        auto real_end = make_real_part_iterator_end(data);

        auto imag_begin = make_imag_part_iterator(data);
        auto imag_end = make_imag_part_iterator_end(data);

        double real_sum = std::accumulate(real_begin, real_end, 0.0);
        double imag_sum = std::accumulate(imag_begin, imag_end, 0.0);

        std::cout << "Real sum: " << real_sum << '\n';
        std::cout << "Imag sum: " << imag_sum << "\n\n";

        FAIL_IF_NOT(real_sum == 10.0);
        FAIL_IF_NOT(imag_sum == 100.0);
    }

    {
        std::cout << "Test 7: empty vector\n";

        std::vector<Complex> empty;

        auto real_begin = make_real_part_iterator(empty);
        auto real_end = make_real_part_iterator_end(empty);

        auto imag_begin = make_imag_part_iterator(empty);
        auto imag_end = make_imag_part_iterator_end(empty);

        FAIL_IF_NOT(real_begin == real_end);
        FAIL_IF_NOT(imag_begin == imag_end);
        FAIL_IF_NOT(real_end - real_begin == 0);
        FAIL_IF_NOT(imag_end - imag_begin == 0);

        std::cout << "Empty vector test passed.\n\n";
    }

    std::cout << "All ComplexNumberIterator tests passed.\n";

    return 0;
}
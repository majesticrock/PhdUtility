#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_INTERNAL_FUNCTIONS_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_INTERNAL_FUNCTIONS_HPP
#include "UnderlyingRealType.hpp"

#include <Eigen/Dense>

#include <charconv>
#include <iostream>

namespace mrock::iEoM {
/**
 * @brief Exception thrown when a matrix contains negative eigenvalues.
 *
 * Carries the most negative eigenvalue found during validation of an
 * eigenvalue vector or diagonal matrix.
 *
 * @tparam RealType Numeric type of the eigenvalue.
 */
template <class RealType>
class MatrixIsNegativeException : public std::runtime_error {
private:
    template <class... Args>
    std::string better_to_string(Args&&... format_args) {
        std::array<char, 32> str;
        auto result = std::to_chars(str.data(), str.data() + str.size(), std::forward<Args>(format_args)...);

        if (result.ec == std::errc())
            return std::string(str.data(), result.ptr - str.data());
        else
            return std::make_error_code(result.ec).message();
    }

public:
    RealType negative_eigenvalue{};

    MatrixIsNegativeException(RealType _negative_eigenvalue, const std::string& name = "M")
        : std::runtime_error("The matrix " + name + " is negative! Most negative eigenvalue = " +
                             better_to_string(_negative_eigenvalue, std::chars_format::scientific, 6)),
          negative_eigenvalue(_negative_eigenvalue){};
};

namespace detail {
enum class iEoM_operation { NONE, INVERSE, SQRT, INVERSE_SQRT };

/**
 * @brief Internal numerical helper for matrix eigenvalue filtering and transformation.
 *
 * Provides precision thresholds, negative-value checking, and
 * eigenvalue operations such as inverse, square root, and inverse square root.
 *
 * @tparam NumberType Scalar type used by the matrices and vectors.
 */
template <class NumberType>
struct iEoM_internal {
    using RealType = UnderlyingRealType_t<NumberType>;
    using RealVector = Eigen::Vector<RealType, Eigen::Dynamic>;

    const RealType _sqrt_precision{1e-6};
    const RealType _precision{1e-12};

    const bool _negative_matrix_is_error{true};

    constexpr iEoM_internal() = default;
    constexpr iEoM_internal(RealType const& sqrt_precision)
        : _sqrt_precision(sqrt_precision), _precision(RealType{2} * sqrt_precision * sqrt_precision){};
    constexpr iEoM_internal(RealType const& sqrt_precision, bool negative_matrix_is_error)
        : _sqrt_precision(sqrt_precision),
          _precision(RealType{2} * sqrt_precision * sqrt_precision),
          _negative_matrix_is_error(negative_matrix_is_error){};

    /**
     * @brief Check whether the vector contains negative eigenvalues.
     *
     * A value is considered negative when it is below the allowed negative
     * precision threshold.
     *
     * @param vector Eigenvalue vector to inspect.
     * @return true if any entry is below -_sqrt_precision.
     */
    inline bool contains_negative(const RealVector& vector) const { return (vector.array() < -_sqrt_precision).any(); };

    /**
     * @brief Apply an eigenvalue operation to a positive semidefinite vector.
     *
     * @tparam option Operation to apply: NONE, INVERSE, SQRT, or INVERSE_SQRT.
     * @tparam pseudo_inverse When true, values below threshold are set to zero.
     * @param evs Vector of eigenvalues to modify in place.
     * @param name Optional name used in diagnostics when a negative eigenvalue is found.
     *
     * This function also verifies that the input vector is non-negative and
     * optionally throws a MatrixIsNegativeException if negative values are detected.
     */
    template <iEoM_operation option, bool pseudo_inverse = true>
    inline void apply_matrix_operation(RealVector& evs, const std::string& name = "M") const {
        if (contains_negative(evs)) {
            if (_negative_matrix_is_error) {
                throw MatrixIsNegativeException<RealType>(evs.minCoeff(), name);
            } else {
                std::cerr << "Warning: The matrix " << name << " is negative with min(ev) = " << evs.minCoeff()
                          << std::endl;
            }
        }

        for (auto& ev : evs) {
            if (ev < _sqrt_precision) {
                if constexpr (pseudo_inverse) {
                    ev = RealType{};
                } else {
                    if constexpr (option == iEoM_operation::INVERSE) {
                        ev = (1. / _precision);
                    } else if constexpr (option == iEoM_operation::SQRT) {
                        ev = _sqrt_precision;
                    } else if constexpr (option == iEoM_operation::INVERSE_SQRT) {
                        ev = (1. / _sqrt_precision);
                    }
                }
            } else {
                if constexpr (option == iEoM_operation::INVERSE) {
                    ev = 1. / ev;
                } else if constexpr (option == iEoM_operation::SQRT) {
                    ev = std::sqrt(ev);
                } else if constexpr (option == iEoM_operation::INVERSE_SQRT) {
                    ev = 1. / std::sqrt(ev);
                }
            }
        }
    };
};
}  // namespace detail
}  // namespace mrock::iEoM
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_INTERNAL_FUNCTIONS_HPP

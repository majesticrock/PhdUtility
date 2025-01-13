#pragma once
#include "../../better_to_string.hpp"
#include "../../UnderlyingFloatingPoint.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace mrock::utility::Numerics::iEoM {
	template<class RealType>
	class MatrixIsNegativeException : public std::runtime_error {
	public:
		RealType negative_eigenvalue{};
		MatrixIsNegativeException(RealType _negative_eigenvalue, const std::string& name = "M")
			: std::runtime_error("The matrix " + name + " is negative! Most negative eigenvalue = "
				+ mrock::utility::better_to_string(_negative_eigenvalue, std::chars_format::scientific, 6)),
			negative_eigenvalue(_negative_eigenvalue)
		{};
	};

	enum ieom_operation { IEOM_NONE, IEOM_INVERSE, IEOM_SQRT, IEOM_INVERSE_SQRT };

	template<class NumberType>
	struct ieom_internal {
		using RealType = UnderlyingFloatingPoint_t<NumberType>;
		using RealVector = Eigen::Vector<RealType, Eigen::Dynamic>;

		const RealType _sqrt_precision{ 1e-6 };
		const RealType _precision{ 1e-12 };

		const bool _negative_matrix_is_error{ true };

		constexpr ieom_internal() = default;
		constexpr ieom_internal(RealType const& sqrt_precision)
			: _sqrt_precision(sqrt_precision), _precision(sqrt_precision* sqrt_precision) {};
		constexpr ieom_internal(RealType const& sqrt_precision, bool negative_matrix_is_error)
			: _sqrt_precision(sqrt_precision), _precision(sqrt_precision* sqrt_precision), _negative_matrix_is_error(negative_matrix_is_error) {};

		template <class EigenMatrixType>
		inline EigenMatrixType remove_noise(EigenMatrixType const& matrix) const {
			return matrix.unaryExpr([this](typename EigenMatrixType::Scalar const& val) {
				return (abs(val) < this->_precision ? typename EigenMatrixType::Scalar{} : val);
				});
		};

		inline bool contains_negative(const RealVector& vector) const {
			return (vector.array() < -_sqrt_precision).any();
		};

		/* Takes a positive semidefinite vector (the idea is that this contains eigenvalues) and applies an operation on it
		* 0: Correct for negative eigenvalues
		* 1: Compute the pseudoinverse
		* 2: Compute the square root
		* 3: Compute the pseudoinverse square root
		*/
		template<ieom_operation option, bool pseudo_inverse = true>
		inline void apply_matrix_operation(RealVector& evs, const std::string& name = "M") const {
			if (contains_negative(evs)) {
				if (_negative_matrix_is_error) {
					throw MatrixIsNegativeException<RealType>(evs.minCoeff(), name);
				}
				else {
					std::cerr << "Warning: The matrix " << name << " is negative with min(ev) = " << evs.minCoeff() << std::endl;
				}
			}

			for (auto& ev : evs)
			{
				if (ev < _sqrt_precision) {
					if constexpr (pseudo_inverse) {
						ev = RealType{};
					}
					else {
						if constexpr (option == IEOM_INVERSE) {
							ev = (1. / _precision);
						}
						else if constexpr (option == IEOM_SQRT) {
							ev = _sqrt_precision;
						}
						else if constexpr (option == IEOM_INVERSE_SQRT) {
							ev = (1. / _sqrt_precision);
						}
					}
				}
				else {
					if constexpr (option == IEOM_INVERSE) {
						ev = 1. / ev;
					}
					else if constexpr (option == IEOM_SQRT) {
						ev = sqrt(ev);
					}
					else if constexpr (option == IEOM_INVERSE_SQRT) {
						ev = 1. / sqrt(ev);
					}
				}
			}
		};
	};
}
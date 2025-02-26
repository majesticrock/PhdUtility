/**
 * @file IterativeSolver.hpp
 * @brief Provides the IterativeSolver template class and related functions for self-consistency computations.
 */

#ifndef _ITERATIVE_SOLVER_HPP
#define _ITERATIVE_SOLVER_HPP

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <limits>
#include <chrono>
#include "../UnderlyingFloatingPoint.hpp"

namespace mrock::utility::Selfconsistency {

	/**
     * @brief Precision constants for different real types.
     */
	template<class RealType> constexpr RealType PRECISION = 1e2 * std::numeric_limits<RealType>::epsilon();
	template<> constexpr float PRECISION<float> = 1e1f * std::numeric_limits<float>::epsilon();

	/**
     * @brief Struct to define debug policies for the iterative solver.
     */
	struct DebugPolicy {
		bool convergenceWarning{ true };
		bool printSteps{ false };
	};
	constexpr DebugPolicy WarnNoConvergence{ true, false };
	constexpr DebugPolicy NoWarning{ false, false };
	constexpr DebugPolicy PrintEverything{ true, true };

	/**
     * @brief Struct to hold convergence information.
     * 
     * @tparam RealType The type of the real number.
     */
	template<class RealType>
	struct ConvergenceInfo {
		RealType error{ 100 };
		bool converged{ false };
		explicit operator bool() const {
			return this->converged;
		};
	};

	/**
     * @brief Overload of the stream insertion operator for ConvergenceInfo.
     * 
     * @tparam RealType The type of the real number.
     * @param os The output stream.
     * @param _info The ConvergenceInfo object.
     * @return The output stream.
     */
	template<class RealType>
	std::ostream& operator<<(std::ostream& os, const ConvergenceInfo<RealType>& _info)
	{
		os << (_info ? "C" : "No c") << "onvergence achieved with error = " << _info.error;
		return os;
	}

	/**
     * @brief Template class for performing iterative self-consistency computations.
     * 
     * @tparam DataType The type of the data.
     * @tparam Model The type of the model.
     * @tparam SelfconsistencyAttributes The type of the self-consistency attributes.
     * @tparam debugPolicy The debug policy to use. Defaults to WarnNoConvergence.
     */
	template <class DataType, class Model, class SelfconsistencyAttributes, const DebugPolicy& debugPolicy = WarnNoConvergence>
	class IterativeSolver {
	protected:
		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;
		using RealType = UnderlyingFloatingPoint_t<DataType>;

		Model* _model{};
		SelfconsistencyAttributes* _attr{};
		const RealType _precision{ PRECISION<RealType> };
		const unsigned int NUMBER_OF_PARAMETERS{};

		/**
         * @brief Checks for sign flipping behavior in the parameter vector.
         * 
         * @param x0 The parameter vector.
         * @return True if sign flipping behavior is detected, false otherwise.
         */
		inline bool has_sign_flipping_behaviour(const ParameterVector& x0) {
			for (unsigned int j = 0U; j < this->NUMBER_OF_PARAMETERS; ++j)
			{
				if (abs(x0[j]) > 1e1 * _precision) {
					if (abs((x0[j] + (*this->_attr)[j]) / x0[j]) < _precision) {
						return true;
					}
				}
			}
			return false;
		}

		/**
         * @brief Performs the iterative procedure for self-consistency.
         * 
         * @param MAX_STEPS The maximum number of steps.
         * @return The convergence information.
         */
		ConvergenceInfo<RealType> procedure_iterative(const unsigned int MAX_STEPS)
		{
			ConvergenceInfo<RealType> convergence;

			ParameterVector f0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
			std::copy(this->_attr->begin(), this->_attr->end(), f0.begin());
			ParameterVector x0{ f0 };

			if (debugPolicy.printSteps) {
				std::cout << "-1:\t" << std::scientific << std::setprecision(4) << "\n" << x0.transpose() << std::endl;
			}

			unsigned int iterNum = 0U;
			while (iterNum < MAX_STEPS && convergence.error > _precision) {
				this->_model->iteration_step(x0, f0);
				if (has_sign_flipping_behaviour(x0)) {
					if constexpr (debugPolicy.convergenceWarning) {
						std::cerr << "Sign flipper for " << this->_model->info() << std::endl;
					}
					return { 2 * x0.norm(), false };
				}
				convergence.error = f0.norm() / x0.norm();
				std::copy(this->_attr->begin(), this->_attr->end(), x0.begin());

				if constexpr (debugPolicy.printSteps) {
					std::cout << iterNum << ":\t" << std::scientific << std::setprecision(4) << "\n" << x0.transpose() << std::endl;
					std::cout << "Error:  " << std::setprecision(12) << convergence.error << "\n";
				}
				++iterNum;
			}
			if constexpr (debugPolicy.printSteps) {
				std::cout << "Finished iterative procedure!" << std::endl;
			}
			if (iterNum >= MAX_STEPS) {
				return convergence;
			}
			convergence.converged = true;
			return convergence;
		}

	public:
		/**
         * @brief Computes the self-consistency attributes.
         * 
         * @param print_time Whether to print the computation time.
         * @param MAX_STEPS The maximum number of steps.
         * @return The self-consistency attributes.
         */
		virtual const SelfconsistencyAttributes& compute(bool print_time = false, const unsigned int MAX_STEPS = 1500)
		{
			std::chrono::time_point begin = std::chrono::steady_clock::now();
			this->_attr->converged = true;
			auto _info = this->procedure_iterative(MAX_STEPS);
			if (!_info) {
				if constexpr (debugPolicy.convergenceWarning) {
					std::cerr << "For " << this->_model->info() << ": " << _info << std::endl;
				}
			}

			if (print_time) {
				std::chrono::time_point end = std::chrono::steady_clock::now();
				std::cout << "Time for self-consistency computations: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			}

			return *this->_attr;
		}

		virtual ~IterativeSolver() = default;
		IterativeSolver() = delete;

		/**
         * @brief Constructs an IterativeSolver with the given model and attributes.
         * 
         * @param model_ptr Pointer to the model.
         * @param attribute_ptr Pointer to the self-consistency attributes.
         */
		IterativeSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr)
			: _model(model_ptr), _attr(attribute_ptr), NUMBER_OF_PARAMETERS(_attr->size()) {}
		
		/**
         * @brief Constructs an IterativeSolver with the given model, attributes, and precision.
         * 
         * @param model_ptr Pointer to the model.
         * @param attribute_ptr Pointer to the self-consistency attributes.
         * @param precision The precision for the iterative procedure.
         */
		IterativeSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, const RealType& precision)
			: _model(model_ptr), _attr(attribute_ptr), _precision(precision), NUMBER_OF_PARAMETERS(_attr->size()) {}
	};

	/**
     * @brief Factory function to create an IterativeSolver with a specified debug policy.
     * 
     * @tparam debugPolicy The debug policy to use.
     * @tparam DataType The type of the data.
     * @tparam Model The type of the model.
     * @tparam SelfconsistencyAttributes The type of the self-consistency attributes.
     * @param model_ptr Pointer to the model.
     * @param attribute_ptr Pointer to the self-consistency attributes.
     * @param precision The precision for the iterative procedure.
     * @return An IterativeSolver object.
     */
	template <const DebugPolicy& debugPolicy, class DataType, class Model, class SelfconsistencyAttributes>
	auto make_iterative(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, const UnderlyingFloatingPoint_t<DataType>& precision = PRECISION<UnderlyingFloatingPoint_t<DataType>>)
	{
		return IterativeSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>(model_ptr, attribute_ptr, precision);
	}

	/**
     * @brief Factory function to create an IterativeSolver with the default debug policy.
     * 
     * @tparam DataType The type of the data.
     * @tparam Model The type of the model.
     * @tparam SelfconsistencyAttributes The type of the self-consistency attributes.
     * @param model_ptr Pointer to the model.
     * @param attribute_ptr Pointer to the self-consistency attributes.
     * @param precision The precision for the iterative procedure.
     * @return An IterativeSolver object.
     */
	template <class DataType, class Model, class SelfconsistencyAttributes>
	auto make_iterative(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, const UnderlyingFloatingPoint_t<DataType>& precision = PRECISION<UnderlyingFloatingPoint_t<DataType>>)
	{
		return IterativeSolver<DataType, Model, SelfconsistencyAttributes, WarnNoConvergence>(model_ptr, attribute_ptr, precision);
	}
}
#endif // _ITERATIVE_SOLVER_HPP
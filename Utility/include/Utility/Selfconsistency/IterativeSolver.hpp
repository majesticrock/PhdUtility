#pragma once
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <limits>
#include <chrono>
#include "../UnderlyingFloatingPoint.hpp"

namespace Utility::Selfconsistency {
	template<class RealType>
	constexpr RealType PRECISION = 1e2 * std::numeric_limits<RealType>::epsilon();
	template<> constexpr float PRECISION<float> = 1e1f * std::numeric_limits<float>::epsilon();

	struct DebugPolicy {
		bool convergenceWarning{ true };
		bool printSteps{ false };
	};
	constexpr DebugPolicy WarnNoConvergence{ true, false };
	constexpr DebugPolicy NoWarning{ false, false };
	constexpr DebugPolicy PrintEverything{ true, true };

	template<class RealType>
	struct ConvergenceInfo {
		RealType error{};
		bool converged{};
		explicit operator bool() const {
			return this->converged;
		};
	};

	template<class RealType>
	std::ostream& operator<<(std::ostream& os, const ConvergenceInfo<RealType>& _info)
	{
		os << (_info ? "C" : "No c") << "onvergence achieved with error = " << _info.error;
		return os;
	}

	template <class DataType, class Model, class SelfconsistencyAttributes, const DebugPolicy& debugPolicy = WarnNoConvergence>
	class IterativeSolver {
	protected:
		Model* _model{};
		SelfconsistencyAttributes* _attr{};
		const size_t NUMBER_OF_PARAMETERS;

		using ParameterVector = Eigen::Vector<DataType, Eigen::Dynamic>;
		using RealType = UnderlyingFloatingPoint_t<DataType>;

		inline bool hasSignFlippingBehaviour(const ParameterVector& x0) {
			for (size_t j = 0U; j < this->NUMBER_OF_PARAMETERS; ++j)
			{
				if (abs(x0[j]) > 1e1 * PRECISION<RealType>) {
					if (abs((x0[j] + (*this->_attr)[j]) / x0[j]) < PRECISION<RealType>) {
						return true;
					}
				}
			}
			return false;
		};

		ConvergenceInfo<RealType> procedureIterative(const size_t MAX_STEPS, RealType precision = PRECISION<RealType>)
		{
			RealType error{ 100 };

			ParameterVector f0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
			std::copy(this->_attr->begin(), this->_attr->end(), f0.begin());
			ParameterVector x0{ f0 };

			if (debugPolicy.printSteps) {
				std::cout << "-1:\t" << std::scientific << std::setprecision(4) << "\n" << x0.transpose() << std::endl;
			}

			size_t iterNum = 0U;
			while (iterNum < MAX_STEPS && error > precision) {
				this->_model->iterationStep(x0, f0);
				if (hasSignFlippingBehaviour(x0)) {
					if constexpr (debugPolicy.convergenceWarning) {
						std::cerr << "Sign flipper for " << this->_model->info() << std::endl;
					}
					return { 2 * x0.norm(), false };
				}
				error = f0.norm();
				std::copy(this->_attr->begin(), this->_attr->end(), x0.begin());

				if constexpr (debugPolicy.printSteps) {
					std::cout << iterNum << ":\t" << std::scientific << std::setprecision(4) << "\n" << x0.transpose() << std::endl;
					std::cout << "Error:  " << error << "\n";
				}
				++iterNum;
			}
			if constexpr (debugPolicy.printSteps) {
				std::cout << "Finished iterative procedure!" << std::endl;
			}
			if (iterNum >= MAX_STEPS) {
				return { error, false };
			}
			return { error, true };
		};
	public:
		virtual const SelfconsistencyAttributes& compute(bool print_time = false, const size_t MAX_STEPS = 1500)
		{
			std::chrono::time_point begin = std::chrono::steady_clock::now();
			this->_attr->converged = true;
			auto _info = this->procedureIterative(MAX_STEPS);
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
		};

		virtual ~IterativeSolver() = default;
		IterativeSolver() = delete;
		IterativeSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr)
			: _model(model_ptr), _attr(attribute_ptr), NUMBER_OF_PARAMETERS(_attr->size()) {};
	};

	template <const DebugPolicy& debugPolicy, class DataType, class Model, class SelfconsistencyAttributes>
	auto make_iterative(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr)
	{
		return IterativeSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>(model_ptr, attribute_ptr);
	}

	template <class DataType, class Model, class SelfconsistencyAttributes>
	auto make_iterative(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr)
	{
		return IterativeSolver<DataType, Model, SelfconsistencyAttributes, WarnNoConvergence>(model_ptr, attribute_ptr);
	}
}
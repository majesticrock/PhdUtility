#ifndef _BROYDEN_SOLVER_HPP
#define _BROYDEN_SOLVER_HPP

#include "IterativeSolver.hpp"
#include "../Numerics/Roots/BroydensMethodEigen.hpp"

namespace mrock::Utility::Selfconsistency {
	template <class DataType, class Model, class SelfconsistencyAttributes, const DebugPolicy& debugPolicy = WarnNoConvergence>
	class BroydenSolver : public IterativeSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>
	{
	private:
		const unsigned int _MaxPreBroydenIterations;
        using _parent = IterativeSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>;
        using ParameterVector = typename _parent::ParameterVector;
        using RealType = typename _parent::RealType;
	public:
		virtual const SelfconsistencyAttributes& compute(bool print_time=false, const unsigned int MAX_STEPS = 400)
		{
            std::chrono::time_point begin = std::chrono::steady_clock::now();
			this->_parent::procedureIterative(_MaxPreBroydenIterations);

			std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
                this->_model->iterationStep(x, F);
			};

			ParameterVector x0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
			std::copy(this->_attr->begin(), this->_attr->end(), x0.begin());
			mrock::Utility::Numerics::Roots::BroydensMethodEigen<DataType, -1> broyden_solver;

			if (!broyden_solver.compute(func, x0, MAX_STEPS)) {
				if (debugPolicy.convergenceWarning) {
					std::cerr << std::fixed << std::setprecision(8) << "No convergence for " << this->_model->info() << std::endl;
				}
			}
			else {
				this->_attr->converged = true;
			}

			if (debugPolicy.printSteps) {
				ParameterVector f0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
				this->_model->iterationStep(x0, f0);
				std::cout << this->_model->info() << "\n";
				std::cout << "x0 = (";
				for (const auto& x : x0)
				{
					std::cout << " " << x << " ";
				}
				std::cout << ")\nf0 = (";
				for (const auto& f : f0)
				{
					std::cout << " " << f << " ";
				}
				std::cout << ")\n -> |f0| = " << std::scientific << std::setprecision(8) << f0.norm() << std::endl;
			}

            if (print_time) {
				std::chrono::time_point end = std::chrono::steady_clock::now();
				std::cout << "Time for self-consistency computations: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
                ParameterVector f0{ ParameterVector::Zero(this->NUMBER_OF_PARAMETERS) };
				this->_model->iterationStep(x0, f0);
				std::cout << "Convergence achieved up to |f0| = " << std::scientific << std::setprecision(8) << f0.norm() << std::endl;
			}

			return *this->_attr;
		};

		BroydenSolver() = delete;
		BroydenSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations)
			: _parent(model_ptr, attribute_ptr, 1e-6), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};
		BroydenSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations, const RealType& preBroydenPrecision)
			: _parent(model_ptr, attribute_ptr, preBroydenPrecision), _MaxPreBroydenIterations(MaxPreBroydenIterations) {};
	};

	template <const DebugPolicy& debugPolicy, class DataType, class Model, class SelfconsistencyAttributes>
	auto make_broyden(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations, const UnderlyingFloatingPoint_t<DataType>& preBroydenPrecision = 1e-6)
	{
		return BroydenSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>(model_ptr, attribute_ptr, MaxPreBroydenIterations, preBroydenPrecision);
	}

	template <class DataType, class Model, class SelfconsistencyAttributes>
	auto make_broyden(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations, const UnderlyingFloatingPoint_t<DataType>& preBroydenPrecision = 1e-6)
	{
		return BroydenSolver<DataType, Model, SelfconsistencyAttributes, WarnNoConvergence>(model_ptr, attribute_ptr, MaxPreBroydenIterations, preBroydenPrecision);
	}
}
#endif
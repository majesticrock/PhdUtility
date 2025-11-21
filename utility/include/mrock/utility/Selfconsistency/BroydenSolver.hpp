/**
 * @file BroydenSolver.hpp
 * @brief Provides the BroydenSolver template class for self-consistency computations using Broyden's method.
 */

#ifndef _BROYDEN_SOLVER_HPP
#define _BROYDEN_SOLVER_HPP

#include "IterativeSolver.hpp"
#include "../Numerics/Roots/BroydensMethodEigen.hpp"

namespace mrock::utility::Selfconsistency {

	/**
     * @brief Template class for performing self-consistency computations using Broyden's method.
     * 
     * @tparam DataType The type of the data.
     * @tparam Model The type of the model.
     * @tparam SelfconsistencyAttributes The type of the self-consistency attributes.
     * @tparam debugPolicy The debug policy to use. Defaults to WarnNoConvergence.
     */
	template <class DataType, class Model, class SelfconsistencyAttributes, const DebugPolicy& debugPolicy = WarnNoConvergence>
	class BroydenSolver : public IterativeSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>
	{
	private:
		const unsigned int _MaxPreBroydenIterations;
        using _parent = IterativeSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>;
        using ParameterVector = typename _parent::ParameterVector;
        using RealType = typename _parent::RealType;

		mrock::utility::Numerics::Roots::BroydensMethodEigen<DataType, -1> broyden_solver;
	public:
		/**
         * @brief Computes the self-consistency attributes using Broyden's method.
         * 
         * @param print_time Whether to print the computation time.
         * @param MAX_STEPS The maximum number of steps.
         * @return The self-consistency attributes.
         */
		virtual const SelfconsistencyAttributes& compute(bool print_time=false, const unsigned int MAX_STEPS = 400, const RealType grace = 1e-10)
		{
            std::chrono::time_point begin = std::chrono::steady_clock::now();
			this->_parent::procedure_iterative(_MaxPreBroydenIterations);

			this->x0.setZero(this->NUMBER_OF_PARAMETERS);
			std::copy(this->_attr->begin(), this->_attr->end(), this->x0.begin());

			std::function<void(const ParameterVector&, ParameterVector&)> func = [&](const ParameterVector& x, ParameterVector& F) {
                this->_model->iteration_step(x, F);
			};
			RealType broyden_convergence = broyden_solver.compute(func, this->x0, MAX_STEPS);
			if (broyden_convergence > grace) {
				this->f0.setZero(this->NUMBER_OF_PARAMETERS);
				this->_model->iteration_step(this->x0, this->f0);

				if (debugPolicy.convergenceWarning) {
					std::cerr << std::fixed << std::setprecision(8) << "No convergence for " << this->_model->info() << std::endl;
					std::cerr << std::scientific << "Final |F(x)| = " << this->f0.norm() << " > " << grace << std::endl;
				}
			}
			else {
				this->_attr->converged = true;
			}

			if (debugPolicy.printSteps) {
				this->f0.setZero(this->NUMBER_OF_PARAMETERS);
				this->_model->iteration_step(this->x0, this->f0);
				std::cout << this->_model->info() << "\n";
				std::cout << "x0 = (";
				for (const auto& x : this->x0)
				{
					std::cout << " " << x << " ";
				}
				std::cout << ")\nf0 = (";
				for (const auto& f : this->f0)
				{
					std::cout << " " << f << " ";
				}
				std::cout << ")\n -> |f0| = " << std::scientific << std::setprecision(8) << this->f0.norm() << std::endl;
			}

            if (print_time) {
				std::chrono::time_point end = std::chrono::steady_clock::now();
				std::cout << "Time for self-consistency computations: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
                this->f0.setZero(this->NUMBER_OF_PARAMETERS);
				this->_model->iteration_step(this->x0, this->f0);
				std::cout << "Convergence achieved up to |f0| = " << std::scientific << std::setprecision(8) << this->f0.norm() << std::endl;
			}

			return *this->_attr;
		}

		virtual void free_memory() override {
			_parent::free_memory();
			broyden_solver.free_memory();
		}

		BroydenSolver() = delete;

		/**
         * @brief Constructs a BroydenSolver with the given model, attributes, and maximum pre-Broyden iterations.
         * 
         * @param model_ptr Pointer to the model.
         * @param attribute_ptr Pointer to the self-consistency attributes.
         * @param MaxPreBroydenIterations The maximum number of pre-Broyden iterations.
         */
		BroydenSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations)
			: _parent(model_ptr, attribute_ptr, 1e-6), _MaxPreBroydenIterations(MaxPreBroydenIterations) {}

		/**
         * @brief Constructs a BroydenSolver with the given model, attributes, maximum pre-Broyden iterations, and precision.
         * 
         * @param model_ptr Pointer to the model.
         * @param attribute_ptr Pointer to the self-consistency attributes.
         * @param MaxPreBroydenIterations The maximum number of pre-Broyden iterations.
         * @param preBroydenPrecision The precision for the pre-Broyden iterative procedure.
         */
		BroydenSolver(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations, const RealType& preBroydenPrecision)
			: _parent(model_ptr, attribute_ptr, preBroydenPrecision), _MaxPreBroydenIterations(MaxPreBroydenIterations) {}
	};

	/**
     * @brief Factory function to create a BroydenSolver with a specified debug policy.
     * 
     * @tparam debugPolicy The debug policy to use.
     * @tparam DataType The type of the data.
     * @tparam Model The type of the model.
     * @tparam SelfconsistencyAttributes The type of the self-consistency attributes.
     * @param model_ptr Pointer to the model.
     * @param attribute_ptr Pointer to the self-consistency attributes.
     * @param MaxPreBroydenIterations The maximum number of pre-Broyden iterations.
     * @param preBroydenPrecision The precision for the pre-Broyden iterative procedure.
     * @return A BroydenSolver object.
     */
	template <const DebugPolicy& debugPolicy, class DataType, class Model, class SelfconsistencyAttributes>
	auto make_broyden(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations, const UnderlyingFloatingPoint_t<DataType>& preBroydenPrecision = 1e-6)
	{
		return BroydenSolver<DataType, Model, SelfconsistencyAttributes, debugPolicy>(model_ptr, attribute_ptr, MaxPreBroydenIterations, preBroydenPrecision);
	}

	/**
     * @brief Factory function to create a BroydenSolver with the default debug policy.
     * 
     * @tparam DataType The type of the data.
     * @tparam Model The type of the model.
     * @tparam SelfconsistencyAttributes The type of the self-consistency attributes.
     * @param model_ptr Pointer to the model.
     * @param attribute_ptr Pointer to the self-consistency attributes.
     * @param MaxPreBroydenIterations The maximum number of pre-Broyden iterations.
     * @param preBroydenPrecision The precision for the pre-Broyden iterative procedure.
     * @return A BroydenSolver object.
     */
	template <class DataType, class Model, class SelfconsistencyAttributes>
	auto make_broyden(Model* model_ptr, SelfconsistencyAttributes* attribute_ptr, unsigned int MaxPreBroydenIterations, const UnderlyingFloatingPoint_t<DataType>& preBroydenPrecision = 1e-6)
	{
		return BroydenSolver<DataType, Model, SelfconsistencyAttributes, WarnNoConvergence>(model_ptr, attribute_ptr, MaxPreBroydenIterations, preBroydenPrecision);
	}
}
#endif // _BROYDEN_SOLVER_HPP
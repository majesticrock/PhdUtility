#pragma once
#include <Eigen/Dense>
#include <vector>
#include <array>
#include "../PivotToBlockStructure.hpp"
#include "_internal_functions.hpp"
#include "../Resolvent.hpp"
#include <chrono>

namespace Utility::Numerics::iEoM {

	template<class Derived, class RealType>
	struct XPResolvent {
	public:
		using Matrix = Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<RealType, Eigen::Dynamic>;

		Matrix K_plus, K_minus, L;
		std::vector<std::array<Vector, 2>> starting_states;
		std::vector<Resolvent<RealType, false>> resolvents;

		// Matrix accessors. Boundary checking is handled by Eigen
		inline const RealType& M(int row, int col) const {
			if(row < _hermitian_size)
				return K_plus(row, col);
			return K_minus(row - _hermitian_size, col - _hermitian_size);
		}
		inline RealType& M(int row, int col) {
			if(row < _hermitian_size)
				return K_plus(row, col);
			return K_minus(row -_hermitian_size, col - _hermitian_size);
		}
		inline const RealType& N(int row, int col) const {
			return L(row, col - _hermitian_size);
		}
		inline RealType& N(int row, int col) {
			return L(row, col - _hermitian_size);
		}

		XPResolvent(Derived* derived_ptr, RealType const& sqrt_precision, bool pivot = true, bool negative_matrix_is_error = true)
			: _internal(sqrt_precision, negative_matrix_is_error), _derived(derived_ptr), _pivot(pivot) { };
		XPResolvent(Derived* derived_ptr, RealType const& sqrt_precision, int hermitian_size, int antihermitian_size,
			bool pivot = true, bool negative_matrix_is_error = true)
			: _internal(sqrt_precision, negative_matrix_is_error), _derived(derived_ptr), 
			_hermitian_size(hermitian_size), _antihermitian_size(antihermitian_size), _pivot(pivot) { };
		virtual ~XPResolvent() = default;

		bool dynamic_matrix_is_negative()
		{
			_derived->fill_M();
			if (_internal.contains_negative(K_minus.diagonal()) || _internal.contains_negative(K_plus.diagonal())) {
				return true;
			}
			K_minus = _internal.removeNoise(K_minus);
			if (not matrix_wrapper<RealType>::is_non_negative(K_minus, _internal._sqrt_precision)) {
				return true;
			}
			K_plus = _internal.removeNoise(K_plus);
			if (not matrix_wrapper<RealType>::is_non_negative(K_plus, _internal._sqrt_precision)) {
				return true;
			}
			return false;
		};

		std::vector<ResolventDataWrapper<RealType>> computeCollectiveModes(unsigned int LANCZOS_ITERATION_NUMBER)
		{
			std::chrono::time_point begin = std::chrono::steady_clock::now();
			_derived->fillMatrices();
			_derived->createStartingStates();

			K_minus = _internal.removeNoise(K_minus);
			K_plus = _internal.removeNoise(K_plus);

			std::chrono::time_point end = std::chrono::steady_clock::now();
			std::cout << "Time for filling of M and N: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			Matrix solver_matrix;
			matrix_wrapper<RealType> k_solutions[2];

			omp_set_nested(2);
			Eigen::initParallel();

#pragma omp parallel sections
			{
#pragma omp section
				{
					std::chrono::time_point begin_in = std::chrono::steady_clock::now();
					if (_pivot) {
						k_solutions[0] = Utility::Numerics::matrix_wrapper<RealType>::pivot_and_solve(K_plus);
					}
					else {
						k_solutions[0] = Utility::Numerics::matrix_wrapper<RealType>::only_solve(K_plus);
					}
					this->_internal.template applyMatrixOperation<IEOM_NONE>(k_solutions[0].eigenvalues, "K_+");
					std::chrono::time_point end_in = std::chrono::steady_clock::now();
					std::cout << "Time for solving K_+: "
						<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

					// free the allocated memory
					K_plus.conservativeResize(0, 0);
				}
#pragma omp section
				{
					std::chrono::time_point begin_in = std::chrono::steady_clock::now();
					if (_pivot) {
						k_solutions[1] = Utility::Numerics::matrix_wrapper<RealType>::pivot_and_solve(K_minus);
					}
					else {
						k_solutions[1] = Utility::Numerics::matrix_wrapper<RealType>::only_solve(K_minus);
					}
					this->_internal.template applyMatrixOperation<IEOM_NONE>(k_solutions[1].eigenvalues, "K_-");
					std::chrono::time_point end_in = std::chrono::steady_clock::now();
					std::cout << "Time for solving K_-: "
						<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

					// free the allocated memory
					K_minus.conservativeResize(0, 0);
				}
			}

			/* plus(minus)_index indicates whether the upper left block is for the
			* Hermitian or the anti-Hermitian operators.
			* The default is that the upper left block contains the Hermtian operators,
			* then plus_index = 0 and minus_index = 1
			*/
			auto compute_solver_matrix = [&](size_t plus_index, size_t minus_index) {
				std::chrono::time_point begin_in = std::chrono::steady_clock::now();
				if (minus_index == 0) L.transposeInPlace();
				solver_matrix.resize(k_solutions[plus_index].eigenvalues.rows(), k_solutions[plus_index].eigenvalues.rows());

				Vector K_EV = k_solutions[minus_index].eigenvalues;
				_internal.template applyMatrixOperation<IEOM_INVERSE>(K_EV, plus_index == 1 ? "K_+" : "K_-");
				Matrix buffer_matrix = L * k_solutions[minus_index].eigenvectors;
				Matrix N_new = buffer_matrix * K_EV.asDiagonal() * buffer_matrix.adjoint();

				std::chrono::time_point end_in = std::chrono::steady_clock::now();
				std::cout << "Time for computing N_new: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
				begin_in = std::chrono::steady_clock::now();

				auto n_solution = _pivot ? Utility::Numerics::matrix_wrapper<RealType>::pivot_and_solve(N_new)
					: Utility::Numerics::matrix_wrapper<RealType>::only_solve(N_new);
				_internal.template applyMatrixOperation<IEOM_INVERSE_SQRT>(n_solution.eigenvalues, plus_index == 1 ? "+: N_new" : "-: N_new");

				// Starting here, N_new = 1/sqrt(N_new)
				// I forego another matrix to save some memory
				N_new = n_solution.eigenvectors * n_solution.eigenvalues.asDiagonal() * n_solution.eigenvectors.adjoint();
				for (auto& starting_state : starting_states) {
					starting_state[plus_index].applyOnTheLeft(N_new * L);
				}

				end_in = std::chrono::steady_clock::now();
				std::cout << "Time for adjusting N_new: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;

				begin_in = std::chrono::steady_clock::now();
				buffer_matrix = N_new * k_solutions[plus_index].eigenvectors;
				solver_matrix = _internal.removeNoise((buffer_matrix * k_solutions[plus_index].eigenvalues.asDiagonal() * buffer_matrix.adjoint()).eval());
				end_in = std::chrono::steady_clock::now();
				std::cout << "Time for computing solver_matrix: "
					<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
				}; // end lambda

			begin = std::chrono::steady_clock::now();

			const int N_RESOLVENT_TYPES = starting_states.size();
			resolvents.reserve(2 * N_RESOLVENT_TYPES);

			for (size_t i = 0U; i < 2U; ++i)
			{
				// It is going to compute the anti-Hermitian first
				compute_solver_matrix(i, 1 - i);
				if(resolvents.size() < 2 * N_RESOLVENT_TYPES){
					for (const auto& starting_state : starting_states) {
						resolvents.push_back(Resolvent<RealType, false>(starting_state[i]));
					}
				} else {
					for(int j = 0; j < N_RESOLVENT_TYPES; ++j) {
						resolvents[N_RESOLVENT_TYPES * i + j].setStartingState(starting_states[j][i]);
					}
				}
#pragma omp parallel for
				for (int j = 0; j < N_RESOLVENT_TYPES; ++j) {
					resolvents[N_RESOLVENT_TYPES * i + j].compute(solver_matrix, LANCZOS_ITERATION_NUMBER);
				}
			}

			end = std::chrono::steady_clock::now();
			std::cout << "Time for resolvents: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			std::vector<ResolventDataWrapper<RealType>> ret;
			ret.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				ret.push_back(re.getData());
			}
			return ret;
		};

	private:
		ieom_internal<RealType> _internal;
		Derived* _derived;
		const int _hermitian_size{};
		const int _antihermitian_size{};
		const bool _pivot{};
	};
}
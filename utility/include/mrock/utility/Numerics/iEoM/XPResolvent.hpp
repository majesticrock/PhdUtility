#pragma once
#include "_internal_functions.hpp"
#include "_xp_internal.hpp"
#include "../PivotToBlockStructure.hpp"
#include "../Resolvent.hpp"
#include "../../ConstexprPower.hpp"
#include "../../OutputConvenience.hpp"
#include <Eigen/Dense>
#include <array>
#include <list>
#include <chrono>

namespace mrock::utility::Numerics::iEoM {
	template<class Derived, class RealType, int n_residuals = 0>
	struct XPResolvent {
	public:
		using Matrix = Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<RealType, Eigen::Dynamic>;

		using phase_it = PhaseIterator<RealType>;
		using amplitude_it = AmplitudeIterator<RealType>;

		Matrix K_plus, K_minus, L;
		std::vector<StartingState<RealType>> starting_states;
		std::vector<Resolvent<Matrix, Vector>> resolvents;

	private:
		template<int CheckHermitian = -1>
		std::array<matrix_wrapper<Matrix>, 2> diagonalize_K_matrices() {
			std::chrono::time_point begin = std::chrono::steady_clock::now();
			_derived->fillMatrices();
			_derived->createStartingStates();

			//K_minus = this->_internal.remove_noise(K_minus);
			//K_plus = this->_internal.remove_noise(K_plus);

			if constexpr (CheckHermitian > 0) {
				if ((K_plus - K_plus.adjoint()).norm() > constexprPower<-CheckHermitian, RealType, RealType>(10.)) {
					throw std::runtime_error("K_plus is not Hermitian!");
				}
				if ((K_minus - K_minus.adjoint()).norm() > constexprPower<-CheckHermitian, RealType, RealType>(10.)) {
					throw std::runtime_error("K_minus is not Hermitian!");
				}
			}

			std::chrono::time_point end = std::chrono::steady_clock::now();
			std::cout << "Time for filling of M and N: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			std::array<matrix_wrapper<Matrix>, 2> k_solutions;

			omp_set_nested(2);
			Eigen::initParallel();

#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel sections
#endif
			{
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp section
#endif
				{
					std::chrono::time_point begin_in = std::chrono::steady_clock::now();
					if (_pivot) {
						k_solutions[0] = matrix_wrapper<Matrix>::pivot_and_solve(K_plus);
					}
					else {
						k_solutions[0] = matrix_wrapper<Matrix>::only_solve(K_plus);
					}
					this->_internal.template apply_matrix_operation<IEOM_NONE>(k_solutions[0].eigenvalues, "K_+");
					std::chrono::time_point end_in = std::chrono::steady_clock::now();
					std::cout << "Time for solving K_+: "
						<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
					// free the allocated memory
					K_plus.resize(0, 0);
				}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp section
#endif
				{
					std::chrono::time_point begin_in = std::chrono::steady_clock::now();
					if (_pivot) {
						k_solutions[1] = matrix_wrapper<Matrix>::pivot_and_solve(K_minus);
					}
					else {
						k_solutions[1] = matrix_wrapper<Matrix>::only_solve(K_minus);
					}
					this->_internal.template apply_matrix_operation<IEOM_NONE>(k_solutions[1].eigenvalues, "K_-");
					std::chrono::time_point end_in = std::chrono::steady_clock::now();
					std::cout << "Time for solving K_-: "
						<< std::chrono::duration_cast<std::chrono::milliseconds>(end_in - begin_in).count() << "[ms]" << std::endl;
					// free the allocated memory
					K_minus.resize(0, 0);
				}
			}
			return k_solutions;
		}

		/* plus(minus)_index indicates whether the upper left block is for the
		* Hermitian or the anti-Hermitian operators.
		* The default is that the upper left block contains the Hermtian operators,
		* then plus_index = 0 and minus_index = 1
		*/
		template <size_t plus_index, size_t minus_index, class StateTransformPolicy>
		void compute_solver_matrix_impl(const std::array<matrix_wrapper<Matrix>, 2>& k_solutions, Matrix& solver_matrix, StateTransformPolicy&& transform) 
		{
			std::chrono::time_point begin = std::chrono::steady_clock::now();
			Vector K_EV = k_solutions[minus_index].eigenvalues;
			_internal.template apply_matrix_operation<IEOM_INVERSE>(K_EV, plus_index == 1 ? "K_+" : "K_-");
			decltype(auto) L_view = [this]() -> decltype(auto) {
			        if constexpr (minus_index == 0)
			            return L.transpose();
			        else
			            return (L);
			    }();
			
			const auto buffer_matrix = L_view * k_solutions[minus_index].eigenvectors;
			Matrix N_new = buffer_matrix * K_EV.asDiagonal() * buffer_matrix.adjoint();

			std::chrono::time_point end = std::chrono::steady_clock::now();
			std::cout << "Time for computing N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			begin = std::chrono::steady_clock::now();	
			{
				auto n_solution = _pivot ? matrix_wrapper<Matrix>::pivot_and_solve(N_new) : matrix_wrapper<Matrix>::only_solve(N_new);
				_internal.template apply_matrix_operation<IEOM_INVERSE_SQRT>(n_solution.eigenvalues, plus_index == 1 ? "+: N_new" : "-: N_new");
				// Starting here, N_new = 1/sqrt(N_new)
				// I forego another matrix to save some memory
				N_new.noalias() = n_solution.eigenvectors * n_solution.eigenvalues.asDiagonal() * n_solution.eigenvectors.adjoint();
				for (auto& starting_state : starting_states) {
					if(starting_state[plus_index].size() > 0U) { 
						transform(starting_state[plus_index], N_new);
					}
				}
			}
			end = std::chrono::steady_clock::now();
			std::cout << "Time for adjusting N_new: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			begin = std::chrono::steady_clock::now();
			N_new *= k_solutions[plus_index].eigenvectors;
			solver_matrix.noalias() = N_new * k_solutions[plus_index].eigenvalues.asDiagonal() * N_new.adjoint();
			//_internal.remove_noise_inplace(solver_matrix);
			end = std::chrono::steady_clock::now();
			std::cout << "Time for computing solver_matrix: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
		}

		/* plus(minus)_index indicates whether the upper left block is for the
		* Hermitian or the anti-Hermitian operators.
		* The default is that the upper left block contains the Hermtian operators,
		* then plus_index = 0 and minus_index = 1
		*/
		template <size_t plus_index, size_t minus_index>
		void compute_solver_matrix(const std::array<matrix_wrapper<Matrix>, 2>& k_solutions, Matrix& solver_matrix) 
		{
			compute_solver_matrix_impl<plus_index, minus_index>(k_solutions, solver_matrix,
				[this](Vector& state, const Matrix& N_new) {
					if constexpr (minus_index == 0) {
						state.applyOnTheLeft(L.transpose());
					}
					else {
						state.applyOnTheLeft(L);
					}
					state.applyOnTheLeft(N_new);
				});
		}

		/* plus(minus)_index indicates whether the upper left block is for the
		* Hermitian or the anti-Hermitian operators.
		* The default is that the upper left block contains the Hermtian operators,
		* then plus_index = 0 and minus_index = 1
		*/
		template <size_t plus_index, size_t minus_index>
		void compute_solver_matrix(const std::array<matrix_wrapper<Matrix>, 2>& k_solutions, Matrix& solver_matrix, Matrix& state_transform) 
		{
			compute_solver_matrix_impl<plus_index, minus_index>(k_solutions, solver_matrix,
				[this, &state_transform](Vector& state, const Matrix& N_new) {
					if constexpr (minus_index == 0) {
						state_transform.noalias() = N_new * L.transpose();
					}
					else {
						state_transform.noalias() = N_new * L;
					}
					state.applyOnTheLeft(state_transform);
				});
		}

	public:
		// Returns a starting state object with an empty phase part and a zero-initialized amplitude part of size /size/
		static StartingState<RealType> OnlyAmplitude(Eigen::Index size, std::string const& name="") {
            return StartingState<RealType>{ Vector{}, Vector::Zero(size), name };
        }
		// Returns a starting state object with an empty amplitude part and a zero-initialized phase part of size /size/
        static StartingState<RealType> OnlyPhase(Eigen::Index size, std::string const& name="") {
            return StartingState<RealType>{ Vector::Zero(size), Vector{}, name };
        }

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
			K_minus = _internal.remove_noise(K_minus);
			if (! matrix_wrapper<Matrix>::is_non_negative(K_minus, _internal._sqrt_precision)) {
				return true;
			}
			K_plus = _internal.remove_noise(K_plus);
			if (! matrix_wrapper<Matrix>::is_non_negative(K_plus, _internal._sqrt_precision)) {
				return true;
			}
			return false;
		};

		template<int CheckHermitian = -1>
		std::vector<resolvent_details::ResolventDataWrapper<RealType>> compute_collective_modes(unsigned int LANCZOS_ITERATION_NUMBER)
		{
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix;

			std::chrono::time_point begin = std::chrono::steady_clock::now();
			const int N_RESOLVENT_TYPES = total_size(starting_states);
			resolvents.resize(N_RESOLVENT_TYPES);
			auto resolvent_it = resolvents.begin();

			compute_solver_matrix<0, 1>(k_solutions, solver_matrix);
			for (phase_it it = phase_it::begin(starting_states); it != phase_it::end(starting_states); ++resolvent_it, ++it) {
				resolvent_it->set_starting_state(it->phase_state);
				if(resolvent_it->data.name.empty()) resolvent_it->data.name = "phase_" + it->name;
			}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
			for (int i = 0; i < phase_size(starting_states); ++i) {
				resolvents[i].compute_with_reorthogonalization(solver_matrix, LANCZOS_ITERATION_NUMBER);
			}

			compute_solver_matrix<1, 0>(k_solutions, solver_matrix);
			for (amplitude_it it = amplitude_it::begin(starting_states); it != amplitude_it::end(starting_states); ++resolvent_it, ++it) {
				resolvent_it->set_starting_state(it->amplitude_state);
				if(resolvent_it->data.name.empty()) resolvent_it->data.name = "amplitude_" + it->name;
			}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
			for (int i = phase_size(starting_states); i < phase_size(starting_states) + amplitude_size(starting_states); ++i) {
				resolvents[i].compute_with_reorthogonalization(solver_matrix, LANCZOS_ITERATION_NUMBER);
			}
			std::chrono::time_point end = std::chrono::steady_clock::now();
			std::cout << "Time for resolvents: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			std::vector<resolvent_details::ResolventDataWrapper<RealType>> ret;
			ret.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				ret.push_back(re.get_data());
			}
			return ret;
		};

		template<int CheckHermitian = -1>
		std::pair<std::vector<resolvent_details::ResolventDataWrapper<RealType>>, std::list<resolvent_details::ResidualInformation<RealType, n_residuals>>> compute_collective_modes_with_residuals(unsigned int LANCZOS_ITERATION_NUMBER)
		{
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix, state_transform;
			std::list<resolvent_details::ResidualInformation<RealType, n_residuals>> residual_infos;

			std::chrono::time_point begin = std::chrono::steady_clock::now();
			const int N_RESOLVENT_TYPES = total_size(starting_states);
			resolvents.resize(N_RESOLVENT_TYPES);
			auto resolvent_it = resolvents.begin();

			compute_solver_matrix<0, 1>(k_solutions, solver_matrix, state_transform);
			for (phase_it it = phase_it::begin(starting_states); it != phase_it::end(starting_states); ++resolvent_it, ++it) {
				resolvent_it->set_starting_state(it->phase_state);
				if(resolvent_it->data.name.empty()) resolvent_it->data.name = "phase_" + it->name;
			}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
			for (int i = 0; i < phase_size(starting_states); ++i) {
				auto residual_info = resolvents[i].template compute_with_residuals<n_residuals>(solver_matrix, LANCZOS_ITERATION_NUMBER);
				for (auto& vj : residual_info.eigenvectors) {
					if (vj.empty()) continue;
					Vector buffer = Vector::Zero(vj.size());
					for (size_t idx = 0U; idx < vj.size(); ++idx) {
						buffer(idx) = vj[idx];
					}
					buffer.applyOnTheLeft(state_transform.adjoint());
					vj = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				for (auto& ev :residual_info.eigenvalues) {
					// The eigenvalue is in z^2
					ev = sqrt(ev);
				}
				residual_infos.emplace_back(std::move(residual_info));
			}

			compute_solver_matrix<1, 0>(k_solutions, solver_matrix, state_transform);
			for (amplitude_it it = amplitude_it::begin(starting_states); it != amplitude_it::end(starting_states); ++resolvent_it, ++it) {
				resolvent_it->set_starting_state(it->amplitude_state);
				if(resolvent_it->data.name.empty()) resolvent_it->data.name = "amplitude_" + it->name;
			}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
			for (int i = phase_size(starting_states); i < phase_size(starting_states) + amplitude_size(starting_states); ++i) {
				auto residual_info = resolvents[i].template compute_with_residuals<n_residuals>(solver_matrix, LANCZOS_ITERATION_NUMBER);
				for (auto& vj : residual_info.eigenvectors) {
					if (vj.empty()) continue;
					Vector buffer = Vector::Zero(vj.size());
					for (size_t idx = 0U; idx < vj.size(); ++idx) {
						buffer(idx) = vj[idx];
					}
					buffer.applyOnTheLeft(state_transform.adjoint());
					vj = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				for (auto& ev :residual_info.eigenvalues) {
					// The eigenvalue is in z^2
					ev = sqrt(ev);
				}
				residual_infos.emplace_back(std::move(residual_info));
			}
			std::chrono::time_point end = std::chrono::steady_clock::now();
			std::cout << "Time for resolvents: "
				<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

			std::vector<resolvent_details::ResolventDataWrapper<RealType>> ret;
			ret.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				ret.push_back(re.get_data());
			}
			return {std::move(ret), std::move(residual_infos)};
		};

		template<int CheckHermitian = -1>
		std::pair<resolvent_details::FullDiagonalizationData<RealType, n_residuals>, resolvent_details::FullDiagonalizationData<RealType, n_residuals>> full_diagonalization() {
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix, state_transform;

			std::chrono::time_point begin = std::chrono::steady_clock::now();
			Eigen::SelfAdjointEigenSolver<Matrix> solver;
			
			resolvent_details::FullDiagonalizationData<RealType, n_residuals> phase_data;
			compute_solver_matrix<0, 1>(k_solutions, solver_matrix, state_transform);
			for (phase_it it = phase_it::begin(starting_states); it != phase_it::end(starting_states); ++it) {
				solver.compute(solver_matrix);

				for (size_t i = 0U; i < n_residuals; ++i){
					Vector buffer = solver.eigenvectors().col(i);
					buffer.applyOnTheLeft(state_transform.adjoint());
					phase_data.first_eigenvectors[i] = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				phase_data.eigenvalues = std::vector<RealType>(solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size());
				phase_data.weights.emplace_back(std::vector<RealType>(solver.eigenvalues().size(), RealType{}));
				Eigen::Map<Vector> weight_map(phase_data.weights.back().data(), solver.eigenvalues().size());
				weight_map = solver.eigenvectors().adjoint() * it->phase_state;
				weight_map = weight_map.array().square();
			}


			resolvent_details::FullDiagonalizationData<RealType, n_residuals> amplitude_data;
			compute_solver_matrix<1, 0>(k_solutions, solver_matrix, state_transform);
			for (amplitude_it it = amplitude_it::begin(starting_states); it != amplitude_it::end(starting_states); ++it) {
				solver.compute(solver_matrix);

				for (size_t i = 0U; i < n_residuals; ++i){
					Vector buffer = solver.eigenvectors().col(i);
					buffer.applyOnTheLeft(state_transform.adjoint());
					amplitude_data.first_eigenvectors[i] = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				amplitude_data.eigenvalues = std::vector<RealType>(solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size());
				amplitude_data.weights.emplace_back(std::vector<RealType>(solver.eigenvalues().size(), RealType{}));
				Eigen::Map<Vector> weight_map(amplitude_data.weights.back().data(), solver.eigenvalues().size());
				weight_map = solver.eigenvectors().adjoint() * it->amplitude_state;
				weight_map = weight_map.array().square();
			}

			return {std::move(phase_data), std::move(amplitude_data)};
		}

	private:
		ieom_internal<RealType> _internal;
		Derived* _derived;
		const int _hermitian_size{};
		const int _antihermitian_size{};
		const bool _pivot{};
	};
}
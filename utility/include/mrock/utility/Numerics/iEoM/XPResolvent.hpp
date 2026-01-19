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
	template<class Derived, class RealType, int n_residuals = 0, bool check_qr = false>
	struct XPResolvent {
	public:
		using Matrix = Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<RealType, Eigen::Dynamic>;

		using phase_it = PhaseIterator<RealType>;
		using amplitude_it = AmplitudeIterator<RealType>;

		using TransformQR = std::conditional_t<!check_qr,
								Eigen::CompleteOrthogonalDecomposition<Eigen::Ref<Matrix>>,
								Eigen::CompleteOrthogonalDecomposition<Matrix>>;

		using FullDiagData = resolvent_details::FullDiagonalizationData<RealType, n_residuals>;
		using ResidualData = resolvent_details::ResidualInformation<RealType, n_residuals>;
		using ResolventReturnData = resolvent_details::ResolventDataWrapper<RealType>;

		Matrix K_plus, K_minus, L;
		std::vector<StartingState<RealType>> starting_states;
		std::vector<Resolvent<Matrix, Vector>> resolvents;

	private:
		std::chrono::time_point<std::chrono::steady_clock> begin;
		std::chrono::time_point<std::chrono::steady_clock> end;

		void set_begin() {
			begin = std::chrono::steady_clock::now();
		}
		void print_duration(const char* message, bool reset=true) {
			end = std::chrono::steady_clock::now();
			std::cout << message << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			if (reset) set_begin();
		}

		template<int CheckHermitian = -1>
		std::array<matrix_wrapper<Matrix>, 2> diagonalize_K_matrices() {
			set_begin();
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

			print_duration("Time for filling of M and N: ");
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
					set_begin();
					if (_pivot) {
						k_solutions[0] = matrix_wrapper<Matrix>::pivot_and_solve(K_plus);
					}
					else {
						k_solutions[0] = matrix_wrapper<Matrix>::only_solve(K_plus);
					}
					this->_internal.template apply_matrix_operation<IEOM_NONE>(k_solutions[0].eigenvalues, "K_+");
					print_duration("Time for solving K_+: ", false);
					// free the allocated memory
					K_plus.resize(0, 0);
				}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp section
#endif
				{
					set_begin();
					if (_pivot) {
						k_solutions[1] = matrix_wrapper<Matrix>::pivot_and_solve(K_minus);
					}
					else {
						k_solutions[1] = matrix_wrapper<Matrix>::only_solve(K_minus);
					}
					this->_internal.template apply_matrix_operation<IEOM_NONE>(k_solutions[1].eigenvalues, "K_-");
					print_duration("Time for solving K_-: ", false);
					// free the allocated memory
					K_minus.resize(0, 0);
				}
			}
			return k_solutions;
		}

		/* plus(minus)_index indicates whether the upper left block is for the
		* Hermitian or the anti-Hermitian operators.
		* The default is that the upper left block contains the Hermitian operators,
		* then plus_index = 0 and minus_index = 1
		*/
		template <size_t plus_index, size_t minus_index, class StateTransformPolicy>
		void compute_solver_matrix_impl(const std::array<matrix_wrapper<Matrix>, 2>& k_solutions, Matrix& solver_matrix, StateTransformPolicy&& transform) 
		{
			set_begin();
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

			print_duration("Time for computing N_new: ");
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
			print_duration("Time for adjusting N_new: ");
			N_new *= k_solutions[plus_index].eigenvectors;
			solver_matrix.noalias() = N_new * k_solutions[plus_index].eigenvalues.asDiagonal() * N_new.adjoint();
			//_internal.remove_noise_inplace(solver_matrix);
			print_duration("Time for computing solver_matrix: ");
		}

		/* plus(minus)_index indicates whether the upper left block is for the
		* Hermitian or the anti-Hermitian operators.
		* The default is that the upper left block contains the Hermitian operators,
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
		* The default is that the upper left block contains the Hermitian operators,
		* then plus_index = 0 and minus_index = 1
		*/
		template <size_t plus_index, size_t minus_index>
		void compute_solver_matrix(const std::array<matrix_wrapper<Matrix>, 2>& k_solutions, Matrix& solver_matrix, Matrix& transform_matrix) 
		{
			compute_solver_matrix_impl<plus_index, minus_index>(k_solutions, solver_matrix,
				[this, &transform_matrix](Vector& state, const Matrix& N_new) {
					if constexpr (minus_index == 0) {
						transform_matrix.noalias() = N_new * L.transpose();
					}
					else {
						transform_matrix.noalias() = N_new * L;
					}
					state.applyOnTheLeft(transform_matrix);
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
		std::vector<ResolventReturnData> compute_collective_modes(unsigned int LANCZOS_ITERATION_NUMBER)
		{
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix;

			set_begin();
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
			print_duration("Time for resolvents: ");

			std::vector<ResolventReturnData> ret;
			ret.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				ret.push_back(re.get_data());
			}
			return ret;
		};

		template<int CheckHermitian = -1>
		std::pair<std::vector<ResolventReturnData>, std::list<ResidualData>> compute_collective_modes_with_residuals(unsigned int LANCZOS_ITERATION_NUMBER)
		{
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix, transform_matrix;

			std::pair<std::vector<ResolventReturnData>, std::list<ResidualData>> return_data;
			std::list<ResidualData>& residual_infos = return_data.second;

			set_begin();
			const int N_RESOLVENT_TYPES = total_size(starting_states);
			resolvents.resize(N_RESOLVENT_TYPES);
			auto resolvent_it = resolvents.begin();

			compute_solver_matrix<0, 1>(k_solutions, solver_matrix, transform_matrix);
			{ TransformQR qr(transform_matrix); // Curly braces to free memory after usage
			print_duration("Time for second QR decomp: ");
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
					Eigen::Map<Vector> vj_eigen(vj.data(), vj.size());
					Vector buffer = qr.solve(vj_eigen);
					if constexpr (check_qr) {
						std::cout << "Error of QR result ||AX - B||=" << (transform_matrix * buffer - vj_eigen).norm() << std::endl;
					}
					vj = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				for (auto& ev :residual_info.eigenvalues) {
					// The eigenvalue is in z^2
					ev = sqrt(ev);
				}
				residual_infos.emplace_back(std::move(residual_info));
			} }

			compute_solver_matrix<1, 0>(k_solutions, solver_matrix, transform_matrix);
			{ TransformQR qr(transform_matrix); // Curly braces to free memory after usage
			print_duration("Time for second QR decomp: ");
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
					Eigen::Map<Vector> vj_eigen(vj.data(), vj.size());
					Vector buffer = qr.solve(vj_eigen);
					if constexpr (check_qr) {
						std::cout << "Error of QR result ||AX - B||=" << (transform_matrix * buffer - vj_eigen).norm() << std::endl;
					}
					vj = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				for (auto& ev :residual_info.eigenvalues) {
					// The eigenvalue is in z^2
					ev = sqrt(ev);
				}
				residual_infos.emplace_back(std::move(residual_info));
			} }

			print_duration("Time for resolvents: ");

			return_data.first.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				return_data.first.push_back(re.get_data());
			}
			return return_data;
		};

		template<int CheckHermitian = -1>
		std::pair<FullDiagData, FullDiagData> full_diagonalization() {
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix, transform_matrix;

			Eigen::SelfAdjointEigenSolver<Matrix> solver;
			std::pair<FullDiagData, FullDiagData> return_data;

			FullDiagData& phase_data = return_data.first;
			compute_solver_matrix<0, 1>(k_solutions, solver_matrix, transform_matrix);

			set_begin();
			solver.compute(solver_matrix);
			print_duration("Time for first ED: ");

			{ TransformQR qr(transform_matrix); // Curly braces to free memory after usage
			print_duration("Time for first QR decomp: ");
			size_t n_zero = 0;
			while (n_zero < solver.eigenvalues().size() && std::abs(solver.eigenvalues()(n_zero)) < _internal._precision) {
				++n_zero;
			}
			const size_t n_non_zero = solver.eigenvalues().size() - n_zero;

			for (phase_it it = phase_it::begin(starting_states); it != phase_it::end(starting_states); ++it) {
				for (size_t i = 0; i < n_residuals; ++i){
					Vector buffer = qr.solve(solver.eigenvectors().col(i + n_zero));
					if constexpr (check_qr) {
						std::cout << "Error of QR result ||AX - B||=" << (transform_matrix * buffer - solver.eigenvectors().col(i + n_zero)).norm() << std::endl;
					}
					phase_data.first_eigenvectors[i] = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				phase_data.eigenvalues = std::vector<RealType>(solver.eigenvalues().data() + n_zero, solver.eigenvalues().data() + solver.eigenvalues().size());
				for (auto& ev : phase_data.eigenvalues) {
					ev = sqrt(ev);
				}
				phase_data.weights.emplace_back(std::vector<RealType>(n_non_zero, RealType{}));
				Eigen::Map<Vector> weight_map(phase_data.weights.back().data(), n_non_zero);
				// phase_state is already transformed; this line computes 
				// sum_j <u_j| transform_matrix | original_phase_state> = sum_j <u_j | N^{-1/2} L | original_phase_state>
				weight_map = solver.eigenvectors().rightCols(n_non_zero).adjoint() * it->phase_state;
				weight_map = weight_map.array().square();
			} }

			FullDiagData& amplitude_data = return_data.second;
			compute_solver_matrix<1, 0>(k_solutions, solver_matrix, transform_matrix);

			set_begin();
			solver.compute(solver_matrix);
			print_duration("Time for second ED: ");

			{ TransformQR qr(transform_matrix); // Curly braces to free memory after usage
			size_t n_zero = 0;
			while (n_zero < solver.eigenvalues().size() && std::abs(solver.eigenvalues()(n_zero)) < _internal._precision) {
				++n_zero;
			}
			const size_t n_non_zero = solver.eigenvalues().size() - n_zero;

			print_duration("Time for second QR decomp: ");
			for (amplitude_it it = amplitude_it::begin(starting_states); it != amplitude_it::end(starting_states); ++it) {
				for (size_t i = 0U; i < n_residuals; ++i){
					Vector buffer = qr.solve(solver.eigenvectors().col(i + n_zero));
					if constexpr (check_qr) {
						std::cout << "Error of QR result ||AX - B||=" << (transform_matrix * buffer - solver.eigenvectors().col(i + n_zero)).norm() << std::endl;
					}
					amplitude_data.first_eigenvectors[i] = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
				}
				amplitude_data.eigenvalues = std::vector<RealType>(solver.eigenvalues().data() + n_zero, solver.eigenvalues().data() + solver.eigenvalues().size());
				for (auto& ev : amplitude_data.eigenvalues) {
					ev = sqrt(ev);
				}
				amplitude_data.weights.emplace_back(std::vector<RealType>(n_non_zero, RealType{}));
				Eigen::Map<Vector> weight_map(amplitude_data.weights.back().data(), n_non_zero);
				// amplitude_state is already transformed; this line computes 
				// sum_j <u_j| transform_matrix | original_amplitude_state> = sum_j <u_j | N^{-1/2} L | original_amplitude_state>
				weight_map = solver.eigenvectors().rightCols(n_non_zero).adjoint() * it->amplitude_state;
				weight_map = weight_map.array().square();
			} }

			return return_data;
		}

	private:
		ieom_internal<RealType> _internal;
		Derived* _derived;
		const int _hermitian_size{};
		const int _antihermitian_size{};
		const bool _pivot{};
	};
}
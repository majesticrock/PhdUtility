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
	/**
	 * @struct XPResolvent
	 * @brief Computes resolvent and eigenvalue data for Hermitian and anti-Hermitian dynamic matrices.
	 * 
	 * This struct performs eigenvalue decomposition and resolvent calculations on K_plus and K_minus
	 * matrices, computing collective modes via Lanczos iteration. It supports both standard resolvent
	 * computation and full diagonalization with residual analysis.
	 * See https://doi.org/10.1103/PhysRevB.109.205153 for the derivation and use of the algorithm.
	 * 
	 * @tparam RealType The floating-point type (e.g., double, float)
	 * @tparam Derived The CRTP-derived class providing matrix filling and starting state creation
	 * @tparam n_residuals Number of residual eigenvectors to retain in full diagonalization
	 * @tparam check_qr If true, validates QR decomposition accuracy via error norms
	 * 
	 * @note Uses OpenMP for parallel eigenvalue solving of K_plus and K_minus when not disabled
	 * @note Maintains timing information for performance profiling via set_begin() and print_duration()
	 * @note Supports optional Hermiticity checking for K_plus and K_minus matrices
	 * 
	 * Public Methods:
	 * - compute_collective_modes(): Computes resolvent functions via Lanczos iteration
	 * - compute_collective_modes_with_residuals(): Extends computation with residual eigenvector information
	 * - full_diagonalization(): Performs complete eigenvalue decomposition retaining first n_residuals eigenvectors
	 * - dynamic_matrix_is_negative(): Checks for negative eigenvalues in diagonal matrices
	 * 
	 * Member Variables:
	 * - K_plus: Hermitian matrix block
	 * - K_minus: Anti-Hermitian matrix block  
	 * - L: Coupling matrix between blocks
	 * - starting_states: Initial phase and amplitude states for resolvent computation
	 * - resolvents: Computed resolvent objects
	 */
	template<class Derived, class RealType, int n_residuals = 0, bool check_qr = false>
	struct XPResolvent {
	public:
		using Matrix = Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<RealType, Eigen::Dynamic>;

		using phase_it = PhaseIterator<RealType>;
		using amplitude_it = AmplitudeIterator<RealType>;

		using const_phase_it = ConstPhaseIterator<RealType>;
		using const_amplitude_it = ConstAmplitudeIterator<RealType>;

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

		/* plus(minus)_index indicates whether the upper left block is for the Hermitian or the anti-Hermitian operators.
		* To compute the Hermitian part set plus_index = 0 and minus_index = 1
		* To compute the anti-Hermitian part set plus_index = 1 and minus_index = 0
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

		/* plus(minus)_index indicates whether the upper left block is for the Hermitian or the anti-Hermitian operators.
		* To compute the Hermitian part set plus_index = 0 and minus_index = 1
		* To compute the anti-Hermitian part set plus_index = 1 and minus_index = 0
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

		/* plus(minus)_index indicates whether the upper left block is for the Hermitian or the anti-Hermitian operators.
		* To compute the Hermitian part set plus_index = 0 and minus_index = 1
		* To compute the anti-Hermitian part set plus_index = 1 and minus_index = 0
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

		/**
		 * @brief Determines the number of zero eigenvalues to exclude from a self-adjoint eigen solver.
		 * 
		 * Counts eigenvalues that fall below the internal precision threshold, treating them as
		 * numerically zero. However, at least one eigenvector from the nullspace is preserved by
		 * decrementing the count if any zero eigenvalues were found.
		 * 
		 * @param solver The Eigen::SelfAdjointEigenSolver containing the computed eigenvalues and eigenvectors.
		 * 
		 * @return The number of zero eigenvalues to exclude, with a minimum of zero ensuring
		 *         at least one nullspace vector is retained when zero eigenvalues exist.
		 * 
		 * @note This function assumes eigenvalues are sorted in ascending order, which Eigen does by default.
		 * @note The precision threshold is obtained from _internal._precision.
		 */
		size_t get_number_of_zero_eigenvalues(const Eigen::SelfAdjointEigenSolver<Matrix>& solver) const noexcept {
			size_t n_zero = 0;
			while (n_zero < solver.eigenvalues().size() && solver.eigenvalues()(n_zero) < _internal._precision) {
				++n_zero;
			}
			if (n_zero > 0U) --n_zero; // Keep one vector from the nullspace
			return n_zero;
		}

		
		/**
		 * @brief Sets full diagonal data for the resolvent calculation.
		 * 
		 * This function prepares and organizes eigenvalue decomposition data along with computed weights for a given state. 
		 * It performs QR-based transformations on eigenvectors and computes weights as squared projections of the state onto
		 * the eigenvector subspace.
		 * TODO: The transformation of the eigenvectors is yet to be published.
		 * 
		 * @param solver A SelfAdjointEigenSolver containing the eigendecomposition of the system matrix.
		 * @param qr A QR decomposition object used to solve linear systems for eigenvector transformation.
		 * @param n_zero The index offset indicating the number of zero/excluded eigenvalues.
		 * @param n_non_zero The number of non-zero eigenvalues and corresponding eigenvectors to process.
		 * @param transform_matrix The transformation matrix used for QR error checking (if enabled).
		 * 
		 * @tparam iterator_type The type of state iterator (const_phase_it or const_amplitude_it) used to traverse starting states.
		 * 
		 * @return FullDiagData A structure containing:
		 *         - first_eigenvectors: Transformed eigenvectors obtained via QR solve
		 *         - eigenvalues: Square root (because the algorithm works in omega^2) of the non-zero eigenvalues 
		 * 								(from index n_zero onward)
		 *         - weights: Squared projections of the state onto the eigenvector subspace,
		 *                    computed as |<u_j | N^{-1/2} L | state>|^2
		 * 
		 * @details If check_qr is true (compile-time constant), the residual error of each
		 *          QR solution is printed for debugging purposes.
		 * 
		 * @note The state vector is assumed to be pre-transformed. The weights represent
		 *       probabilities (squared amplitudes) in the eigenvector basis.
		 */
		template<ConstStateIterator iterator_type>
		FullDiagData set_full_diag_data(const Eigen::SelfAdjointEigenSolver<Matrix>& solver, const TransformQR& qr, 
			const size_t& n_zero, const size_t& n_non_zero, const Matrix& transform_matrix) const
		{
			FullDiagData data;
			for (size_t i = 0U; i < n_residuals; ++i){
				Vector buffer = qr.solve(solver.eigenvectors().col(i + n_zero));
				if constexpr (check_qr) {
					std::cout << "Error of QR result ||AX - B||=" << (transform_matrix * buffer - solver.eigenvectors().col(i + n_zero)).norm() << std::endl;
				}
				data.first_eigenvectors[i] = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
			}
			data.eigenvalues = std::vector<RealType>(solver.eigenvalues().data() + n_zero, solver.eigenvalues().data() + solver.eigenvalues().size());
			for (auto& ev : data.eigenvalues) {
				ev = sqrt(ev);
			}
			for (iterator_type it = iterator_type::begin(starting_states); it != iterator_type::end(starting_states); ++it) {
				data.weights.emplace_back(std::vector<RealType>(n_non_zero, RealType{}));
				Eigen::Map<Vector> weight_map(data.weights.back().data(), n_non_zero);
				// state is already transformed; this line computes
				// sum_j <u_j| transform_matrix | original_amplitude_state> = sum_j <u_j | N^{-1/2} L | original_amplitude_state>
				// The analogue also applies to phase states:  sum_j <u_j | N^{-1/2} L^+ | original_phase_state>
				weight_map = solver.eigenvectors().rightCols(n_non_zero).adjoint() * it.state();
				weight_map = weight_map.array().square();
			}
			return data;
		}

		/**
		 * @brief Sets residual data by computing eigenvalues and eigenvectors with the residuals of the Lanczos algorithm.
		 * 			This part is handled by the Resolvent class.
		 * 
		 * This function computes the eigenvalues and eigenvectors of the solver matrix using
		 * Lanczos iteration with residual information. It then transforms the eigenvectors
		 * using QR decomposition and takes the square root of the eigenvalues.
		 * 
		 * @param resolvent Reference to the Resolvent object used to compute eigendecomposition.
		 * @param n_lanczos_iterations Number of Lanczos iterations to perform.
		 * @param qr The QR transformation used to solve for the transformed eigenvectors.
		 * @param transform_matrix The transformation matrix used for QR error checking (if enabled).
		 * 
		 * @return ResidualData Structure containing the transformed eigenvectors and computed eigenvalues.
		 *         - eigenvectors: QR-transformed eigenvectors (empty vectors are skipped)
		 *         - eigenvalues: Square roots of the computed eigenvalues
		 * 
		 * @note If check_qr is true at compile-time, outputs the residual error ||AX - B|| for each eigenvector.
		 * @note Eigenvalues are assumed to be in z^2 form and are converted to z via sqrt().
		 */
		ResidualData set_residual_data(Resolvent<Matrix, Vector>& resolvent, int n_lanczos_iterations, 
			const Matrix& solver_matrix, const TransformQR& qr, const Matrix& transform_matrix) const 
		{
			ResidualData residual_info = resolvent.template compute_with_residuals<n_residuals>(solver_matrix, n_lanczos_iterations);
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
			return residual_info;
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

		// Matrix accessors. Boundary checking is handled by Eigen - if NBEDUG is not defined.
		inline const RealType& M(int row, int col) const {
			if(row < _hermitian_size)
				return K_plus.coeffRef(row, col);
			return K_minus.coeffRef(row - _hermitian_size, col - _hermitian_size);
		}
		inline RealType& M(int row, int col) {
			if(row < _hermitian_size)
				return K_plus.coeffRef(row, col);
			return K_minus.coeffRef(row -_hermitian_size, col - _hermitian_size);
		}
		inline const RealType& N(int row, int col) const {
			return L.coeffRef(row, col - _hermitian_size);
		}
		inline RealType& N(int row, int col) {
			return L.coeffRef(row, col - _hermitian_size);
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
		std::vector<ResolventReturnData> compute_collective_modes(unsigned int n_lanczos_iterations)
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
				resolvents[i].compute_with_reorthogonalization(solver_matrix, n_lanczos_iterations);
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
				resolvents[i].compute_with_reorthogonalization(solver_matrix, n_lanczos_iterations);
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
		std::pair<std::vector<ResolventReturnData>, std::list<ResidualData>> compute_collective_modes_with_residuals(unsigned int n_lanczos_iterations)
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
			{ 
				TransformQR qr(transform_matrix); // Curly braces to free memory after usage
				print_duration("Time for second QR decomp: ");
				for (phase_it it = phase_it::begin(starting_states); it != phase_it::end(starting_states); ++resolvent_it, ++it) {
					resolvent_it->set_starting_state(it->phase_state);
					if(resolvent_it->data.name.empty()) resolvent_it->data.name = "phase_" + it->name;
				}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
				for (int i = 0; i < phase_size(starting_states); ++i) {
					residual_infos.emplace_back(set_residual_data(resolvents[i], n_lanczos_iterations, solver_matrix, qr, transform_matrix));
				} 
			}

			compute_solver_matrix<1, 0>(k_solutions, solver_matrix, transform_matrix);
			{ 
				TransformQR qr(transform_matrix); // Curly braces to free memory after usage
				print_duration("Time for second QR decomp: ");
				for (amplitude_it it = amplitude_it::begin(starting_states); it != amplitude_it::end(starting_states); ++resolvent_it, ++it) {
					resolvent_it->set_starting_state(it->amplitude_state);
					if(resolvent_it->data.name.empty()) resolvent_it->data.name = "amplitude_" + it->name;
				}
#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#pragma omp parallel for
#endif
				for (int i = phase_size(starting_states); i < phase_size(starting_states) + amplitude_size(starting_states); ++i) {
					residual_infos.emplace_back(set_residual_data(resolvents[i], n_lanczos_iterations, solver_matrix, qr, transform_matrix));
				} 
			}
			print_duration("Time for resolvents: ");

			return_data.first.reserve(resolvents.size());
			for (const auto& re : resolvents)
			{
				return_data.first.push_back(re.get_data());
			}
			return return_data;
		};

		/*
		* Performs a full diagonalization of the dynamic matrices.
		* Saves the first <n_residuals> eigenvectors, all eigenvalues and the weights for each starting state.
		* Returns a pair of FullDiagData objects, the first for the phase states, the second for the amplitude states.
		*/
		template<int CheckHermitian = -1>
		std::pair<FullDiagData, FullDiagData> full_diagonalization() {
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix, transform_matrix;

			Eigen::SelfAdjointEigenSolver<Matrix> solver;
			std::pair<FullDiagData, FullDiagData> return_data;

			compute_solver_matrix<0, 1>(k_solutions, solver_matrix, transform_matrix);
			set_begin();
			solver.compute(solver_matrix);
			print_duration("Time for first ED: ");
			{ // Curly braces to free memory after usage
				TransformQR qr(transform_matrix); 
				print_duration("Time for first QR decomp: ");
				const size_t n_zero = get_number_of_zero_eigenvalues(solver);
				const size_t n_non_zero = solver.eigenvalues().size() - n_zero;
				return_data.first = set_full_diag_data<const_phase_it>(solver, qr, n_zero, n_non_zero, transform_matrix);
				if constexpr (check_qr) {
					for (size_t i = 0U; i < n_residuals; ++i) {
						Eigen::Map<Vector> _eigen(return_data.first.first_eigenvectors[i].data(), return_data.first.first_eigenvectors[i].size());
						const auto reconstructed = transform_matrix * _eigen;
						std::cout << "|solve_matrix * reconstructed eigenvector - eigenvalue * reconstructed eigenvector| = " 
							<< (solver_matrix * reconstructed - return_data.first.eigenvalues[i] * return_data.first.eigenvalues[i] * reconstructed).norm() << std::endl;
					}
				}
			}

			compute_solver_matrix<1, 0>(k_solutions, solver_matrix, transform_matrix);
			set_begin();
			solver.compute(solver_matrix);
			print_duration("Time for second ED: ");
			{ // Curly braces to free memory after usage
				TransformQR qr(transform_matrix); 
				const size_t n_zero = get_number_of_zero_eigenvalues(solver);
				const size_t n_non_zero = solver.eigenvalues().size() - n_zero;

				print_duration("Time for second QR decomp: ");
				return_data.second = set_full_diag_data<const_amplitude_it>(solver, qr, n_zero, n_non_zero, transform_matrix);
				if constexpr (check_qr) {
					for (size_t i = 0U; i < n_residuals; ++i) {
						Eigen::Map<Vector> _eigen(return_data.second.first_eigenvectors[i].data(), return_data.second.first_eigenvectors[i].size());
						const auto reconstructed = transform_matrix * _eigen;
						std::cout << "|solve_matrix * reconstructed eigenvector - eigenvalue * reconstructed eigenvector| = " 
							<< (solver_matrix * reconstructed - return_data.second.eigenvalues[i] * return_data.second.eigenvalues[i] * reconstructed).norm() << std::endl;
					}
				}
			}

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
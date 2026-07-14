#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_XPRESOLVENT_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_XPRESOLVENT_HPP

#ifndef _OPENMP
#define MROCK_IEOM_DO_NOT_PARALLELIZE
#endif

#ifndef MROCK_IEOM_DO_NOT_PARALLELIZE
#include <omp.h>
#endif

#include <Eigen/Dense>
#include <array>
#include <list>
#include <chrono>

#include "detail/internal_functions.hpp"
#include "detail/xp_internal.hpp"
#include "detail/PivotToBlockStructure.hpp"
#include "detail/constexpr_power.hpp"
#include "XPStartingState.hpp"
#include "Resolvent.hpp"

namespace mrock::iEoM {
	/**
	 * @struct XPResolvent
	 * @brief Computes resolvent and eigenvalue data for Hermitian and anti-Hermitian dynamical matrices.
	 * 
	 * This struct performs eigenvalue decomposition and resolvent calculations on K_plus and K_minus
	 * matrices, computing collective modes via Lanczos iteration. It supports both standard resolvent
	 * computation and full diagonalization with residual analysis.
	 * See https://doi.org/10.1103/PhysRevB.109.205153 for the derivation and use of the algorithm.
	 * The full_diagonalization() and compute_collective_modes_with_residuals() methods are specially
	 * designed to identify the operators that excite a given mode ("eigenoperators").
	 * The algorithm is derived in https://doi.org/10.48550/arXiv.2605.20059
	 * 
	 * @tparam RealType The floating-point type (e.g., double, float)
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
	 * - K_plus: Hermitian block of the dynamical matrix
	 * - K_minus: Anti-Hermitian block of the dynamical matrix
	 * - L: norm matrix; couples between blocks
	 * - starting_states: Initial phase and amplitude states for resolvent computation
	 */
	template<class RealType, int n_residuals = 0, bool check_qr = false>
	struct XPResolvent {
	public:
		using ResolventReturnData = ResolventDataWrapper<RealType>;
		using FullDiagData = FullDiagonalizationData<RealType, n_residuals>;
		using ResidualData = ResidualInformation<RealType, n_residuals>;

	protected:
		using Matrix = Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic>;
		using Vector = Eigen::Vector<RealType, Eigen::Dynamic>;

		using phase_it = detail::PhaseIterator<RealType>;
		using amplitude_it = detail::AmplitudeIterator<RealType>;

		using const_phase_it = detail::ConstPhaseIterator<RealType>;
		using const_amplitude_it = detail::ConstAmplitudeIterator<RealType>;

		using StartingState = XPStartingState<RealType>;

		using TransformQR = std::conditional_t<!check_qr,
								Eigen::CompleteOrthogonalDecomposition<Eigen::Ref<Matrix>>,
								Eigen::CompleteOrthogonalDecomposition<Matrix>>;

		Matrix K_plus, K_minus, L;
		std::vector<StartingState> starting_states;
	
	private:
		detail::iEoM_internal<RealType> _internal;
		const int _hermitian_size{};
		const int _antihermitian_size{};
		const bool _pivot{};
		bool compute_transform{true};

		std::chrono::time_point<std::chrono::steady_clock> begin;
		std::chrono::time_point<std::chrono::steady_clock> end;

	protected:
		// The following virtual functions must be implemented by the user
		// If one knows that certain functions are not needed, implementing an empty
		// function suffices.
		// An example would be: fill_M() is only needed by dynamic_matrix_is_negative()
		// Thus, if this function is not called, fill_M() is not required either.

		/**
		 * @brief Must be implemented by the user. Must fill the starting states for the Laczos algorithm
		 */
		virtual void create_starting_states() = 0;
		/**
		 * @brief Must be implemented by the user. Must fill the matrices K_+, K_-, and L
		 */
		virtual void fill_matrices() = 0;
		/**
		 * @brief Must be implemented by the user. Must fill the matrices K_+ and K_-
		 */
		virtual void fill_M() = 0;

		/**
		 * @brief Access an element of the dynamical matrix M.
		 *
		 * @param row Row index in the combined Hermitian/anti-Hermitian matrix.
		 * @param col Column index in the combined Hermitian/anti-Hermitian matrix.
		 * @return Reference to the requested matrix element.
		 *
		 * Boundary checking is handled by Eigen unless NBEDUG is defined.
		 */
		inline const RealType& M(int row, int col) const {
			if(row < _hermitian_size)
				return K_plus.coeffRef(row, col);
			return K_minus.coeffRef(row - _hermitian_size, col - _hermitian_size);
		}
		/**
		 * @brief Access a mutable element of the dynamical matrix M.
		 *
		 * @param row Row index in the combined Hermitian/anti-Hermitian matrix.
		 * @param col Column index in the combined Hermitian/anti-Hermitian matrix.
		 * @return Mutable reference to the requested matrix element.
		 *
		 * Boundary checking is handled by Eigen unless NBEDUG is defined.
		 */
		inline RealType& M(int row, int col) {
			if(row < _hermitian_size)
				return K_plus.coeffRef(row, col);
			return K_minus.coeffRef(row -_hermitian_size, col - _hermitian_size);
		}
		/**
		 * @brief Access an element of the norm matrix N.
		 *
		 * @param row Row index in the norm matrix.
		 * @param col Column index adjusted for the Hermitian block offset.
		 * @return Reference to the requested norm matrix element.
		 */
		inline const RealType& N(int row, int col) const {
			return L.coeffRef(row, col - _hermitian_size);
		}
		/**
		 * @brief Access a mutable element of the norm matrix N.
		 *
		 * @param row Row index in the norm matrix.
		 * @param col Column index adjusted for the Hermitian block offset.
		 * @return Mutable reference to the requested norm matrix element.
		 */
		inline RealType& N(int row, int col) {
			return L.coeffRef(row, col - _hermitian_size);
		}

	private:
		/**
		 * @brief Start or restart the internal profiling timer.
		 *
		 * Records the current steady-clock time into the internal `begin` member.
		 * This is used in conjunction with `print_duration()` to measure elapsed
		 * time for sections of the algorithm for lightweight profiling.
		 */
		void set_begin() {
			begin = std::chrono::steady_clock::now();
		}
		
		/**
		 * @brief Print elapsed time since last call to `set_begin()`.
		 *
		 * @param message A C-string prefix printed before the elapsed time (in milliseconds).
		 * @param reset If true, restarts the timer after printing (default: true).
		 *
		 * This helper is used for lightweight profiling and debugging. Time is measured
		 * using `std::chrono::steady_clock` and printed to `std::cout` in milliseconds.
		 */
		void print_duration(const char* message, bool reset=true) {
			end = std::chrono::steady_clock::now();
			std::cout << message << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
			if (reset) set_begin();
		}

		/**
		 * @brief Diagonalize the K_plus and K_minus dynamical matrices.
		 *
		 * This function delegates to the derived class to fill the matrix blocks and
		 * to create starting states, optionally performs a Hermiticity check (when
		 * `CheckHermitian>0`), and then computes eigen-decompositions of both
		 * `K_plus` and `K_minus`. Eigen-solves are executed in parallel where
		 * permitted. After decomposition the on-device `K_plus` and `K_minus` storage
		 * is released (resized to 0).
		 *
		 * @tparam CheckHermitian If positive, enable an additional Hermiticity
		 *                        verification using a precision threshold of
		 *                        10^-CheckHermitian.
		 * @return An array with two `detail::matrix_wrapper<Matrix>` entries holding
		 *         the eigenvalues and eigenvectors for `K_plus` (index 0) and
		 *         `K_minus` (index 1) respectively.
		 */
		template<int CheckHermitian = -1>
		std::array<detail::matrix_wrapper<Matrix>, 2> diagonalize_K_matrices() {
			set_begin();
			fill_matrices();
			create_starting_states();

			if constexpr (CheckHermitian > 0) {
				if ((K_plus - K_plus.adjoint()).norm() > detail::constexpr_power<-CheckHermitian, RealType, RealType>(10.)) {
					throw std::runtime_error("K_plus is not Hermitian!");
				}
				if ((K_minus - K_minus.adjoint()).norm() > detail::constexpr_power<-CheckHermitian, RealType, RealType>(10.)) {
					throw std::runtime_error("K_minus is not Hermitian!");
				}
			}

			print_duration("Time for filling of M and N: ");
			std::array<detail::matrix_wrapper<Matrix>, 2> k_solutions;

#ifdef _OPENMP
			omp_set_nested(2);
			Eigen::initParallel();
#endif

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
						k_solutions[0] = detail::matrix_wrapper<Matrix>::pivot_and_solve(K_plus);
					}
					else {
						k_solutions[0] = detail::matrix_wrapper<Matrix>::only_solve(K_plus);
					}
					this->_internal.template apply_matrix_operation<detail::iEoM_operation::NONE>(k_solutions[0].eigenvalues, "K_+");
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
						k_solutions[1] = detail::matrix_wrapper<Matrix>::pivot_and_solve(K_minus);
					}
					else {
						k_solutions[1] = detail::matrix_wrapper<Matrix>::only_solve(K_minus);
					}
					this->_internal.template apply_matrix_operation<detail::iEoM_operation::NONE>(k_solutions[1].eigenvalues, "K_-");
					print_duration("Time for solving K_-: ", false);
					// free the allocated memory
					K_minus.resize(0, 0);
				}
			}
			return k_solutions;
		}

		/**
		 * @brief Computes the solver matrix using eigen-decomposed K matrices.
		 *
		 * @tparam plus_index Index of the K matrix corresponding to the Hermitian block in the assembled solver.
		 * @tparam minus_index Index of the K matrix corresponding to the anti-Hermitian block in the assembled solver.
		 * @tparam StateTransformPolicy Transformation functor applied to each starting state.
		 *
		 * @param k_solutions Eigen-decomposed K_plus and K_minus wrapped matrices.
		 * @param solver_matrix Output matrix that will contain the computed solver matrix.
		 * @param transform Functor used to transform each starting state with the computed N_new matrix.
		 *
		 * To compute the Hermitian solver part use plus_index = 0 and minus_index = 1.
		 * To compute the anti-Hermitian solver part use plus_index = 1 and minus_index = 0.
		 */
		template <std::size_t plus_index, std::size_t minus_index, class StateTransformPolicy>
		void compute_solver_matrix_impl(const std::array<detail::matrix_wrapper<Matrix>, 2>& k_solutions, 
			Matrix& solver_matrix, StateTransformPolicy&& transform) 
		{
			set_begin();
			Vector K_EV = k_solutions[minus_index].eigenvalues;
			_internal.template apply_matrix_operation<detail::iEoM_operation::INVERSE>(K_EV, plus_index == 1 ? "K_+" : "K_-");
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
				auto n_solution = _pivot ? detail::matrix_wrapper<Matrix>::pivot_and_solve(N_new) : detail::matrix_wrapper<Matrix>::only_solve(N_new);
				_internal.template apply_matrix_operation<detail::iEoM_operation::INVERSE_SQRT>(n_solution.eigenvalues, plus_index == 1 ? "+: N_new" : "-: N_new");
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
			
			print_duration("Time for computing solver_matrix: ");
		}

		/**
		 * @brief Computes solver_matrix and applies the standard state transformation.
		 *
		 * @tparam plus_index Index of the Hermitian block within the assembled solver.
		 * @tparam minus_index Index of the anti-Hermitian block within the assembled solver.
		 *
		 * @param k_solutions Eigen-decomposed K_plus and K_minus wrapped matrices.
		 * @param solver_matrix Output matrix populated by the solver computation.
		 *
		 * Uses the standard state transform path for phase and amplitude states.
		 */
		template <std::size_t plus_index, std::size_t minus_index>
		void compute_solver_matrix(const std::array<detail::matrix_wrapper<Matrix>, 2>& k_solutions, Matrix& solver_matrix) 
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

		/**
		 * @brief Computes solver_matrix and applies a custom transformation matrix.
		 *
		 * @tparam plus_index Index of the Hermitian block within the assembled solver.
		 * @tparam minus_index Index of the anti-Hermitian block within the assembled solver.
		 *
		 * @param k_solutions Eigen-decomposed K_plus and K_minus wrapped matrices.
		 * @param solver_matrix Output matrix populated by the solver computation.
		 * @param transform_matrix Matrix used to transform each state before applying it to the solver.
		 *
		 * Uses a precomputed transformation matrix to update the solver state vectors.
		 */
		template <std::size_t plus_index, std::size_t minus_index>
		void compute_solver_matrix(const std::array<detail::matrix_wrapper<Matrix>, 2>& k_solutions, Matrix& solver_matrix, Matrix& transform_matrix) 
		{
			compute_transform = true;
			compute_solver_matrix_impl<plus_index, minus_index>(k_solutions, solver_matrix,
				[this, &transform_matrix](Vector& state, const Matrix& N_new) {
					if (!compute_transform) {
						state.applyOnTheLeft(transform_matrix);
						return;
					}

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
		std::size_t get_number_of_zero_eigenvalues(const Eigen::SelfAdjointEigenSolver<Matrix>& solver) const noexcept {
			std::size_t n_zero{};
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
		template<detail::ConstStateIterator iterator_type>
		FullDiagData set_full_diag_data(const Eigen::SelfAdjointEigenSolver<Matrix>& solver, const TransformQR& qr, 
			const std::size_t& n_zero, const std::size_t& n_non_zero, const Matrix& transform_matrix) const
		{
			FullDiagData data;
			for (std::size_t i = 0U; i < n_residuals; ++i){
				Vector buffer = qr.solve(solver.eigenvectors().col(i + n_zero));
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
				// The analog also applies to phase states:  sum_j <u_j | N^{-1/2} L^+ | original_phase_state>
				weight_map = solver.eigenvectors().rightCols(n_non_zero).adjoint() * it.state();
				weight_map = weight_map.array().square();

				const Vector nullspace_weights = solver.eigenvectors().leftCols(n_zero).adjoint() * it.state();
				for (const auto& weight : nullspace_weights) {
					data.weights.back()[0] += weight;
				}
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
				vj = std::vector<RealType>(buffer.data(), buffer.data() + buffer.size());
			}
			for (auto& ev :residual_info.eigenvalues) {
				// The eigenvalue is in z^2
				ev = sqrt(ev);
			}
			return residual_info;
		}

	public:
		/**
		 * @brief Create a starting state containing only an amplitude component.
		 *
		 * @param size Number of amplitude entries to initialize.
		 * @param name Optional name for the starting state.
		 * @return StartingState with an empty phase_state and a zero-initialized amplitude_state.
		 */
		static StartingState OnlyAmplitude(Eigen::Index size, std::string const& name="") {
            return StartingState{ Vector{}, Vector::Zero(size), name };
        }
		/**
		 * @brief Create a starting state containing only a phase component.
		 *
		 * @param size Number of phase entries to initialize.
		 * @param name Optional name for the starting state.
		 * @return StartingState with a zero-initialized phase_state and an empty amplitude_state.
		 */
        static StartingState OnlyPhase(Eigen::Index size, std::string const& name="") {
            return StartingState{ Vector::Zero(size), Vector{}, name };
        }

		/**
		 * @brief Construct an XPResolvent
		 *
		 * @param sqrt_precision Square-root precision threshold for internal checks.
		 * @param pivot Enable matrix pivoting during solves (default: true).
		 * @param negative_matrix_is_error If true, negative diagonal entries are treated as errors.
		 */
		XPResolvent(RealType const& sqrt_precision, bool pivot = true, bool negative_matrix_is_error = true)
			: _internal(sqrt_precision, negative_matrix_is_error), _pivot(pivot) { };

		/**
		 * @brief Construct an XPResolvent with explicit block sizes.
		 *
		 * This constructor allows directly specifying the sizes of the Hermitian and
		 * anti-Hermitian blocks when those are known ahead of time.
		 *
		 * @param sqrt_precision Square-root precision threshold for internal checks.
		 * @param hermitian_size Number of rows/cols in the Hermitian block.
		 * @param antihermitian_size Number of rows/cols in the anti-Hermitian block.
		 * @param pivot Enable matrix pivoting during solves (default: true).
		 * @param negative_matrix_is_error If true, negative diagonal entries are treated as errors.
		 */
		XPResolvent(RealType const& sqrt_precision, int hermitian_size, int antihermitian_size,
			bool pivot = true, bool negative_matrix_is_error = true)
			: _internal(sqrt_precision, negative_matrix_is_error), 
			_hermitian_size(hermitian_size), _antihermitian_size(antihermitian_size), _pivot(pivot) { };

		/**
		 * @brief Virtual default destructor.
		 */
		virtual ~XPResolvent() = default;

		/**
		 * @brief Checks whether the assembled dynamical matrix contains negative eigenvalues,
		 * or its diagonal contains negative numbers.
		 * These must be >= 0 in thermal equilibrium. Thus, if this function returns true,
		 * one can deduce that the system under study is not in thermal equilibrium.
		 * Note that the converse is not true.
		 *
		 * @return true when either K_minus or K_plus contains negative eigenvalues.
		 */
		bool dynamic_matrix_is_negative()
		{
			fill_M();
			if (_internal.contains_negative(K_minus.diagonal()) || _internal.contains_negative(K_plus.diagonal())) {
				return true;
			}
			if (! detail::matrix_wrapper<Matrix>::is_non_negative(K_minus, _internal._sqrt_precision)) {
				return true;
			}
			if (! detail::matrix_wrapper<Matrix>::is_non_negative(K_plus, _internal._sqrt_precision)) {
				return true;
			}
			return false;
		};

		/**
		 * @brief Compute collective-mode resolvents using Lanczos iteration.
		 *
		 * This routine diagonalizes the K_plus and K_minus blocks (optionally checking
		 * Hermiticity controlled by the `CheckHermitian` template parameter), assembles
		 * solver matrices for phase and amplitude channels, then runs Lanczos iterations
		 * to obtain resolvent data for all configured starting states.
		 *
		 * @tparam CheckHermitian If >0, enables runtime Hermiticity checks against precision 10^-CheckHermitian.
		 * @param n_lanczos_iterations Number of Lanczos iterations used for each resolvent.
		 * @return A vector of `ResolventReturnData` containing resolvent results for
		 *         phase followed by amplitude starting states.
		 */
		template<int CheckHermitian = -1>
		std::vector<ResolventReturnData> compute_collective_modes(unsigned int n_lanczos_iterations)
		{
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix;

			set_begin();
			const int N_RESOLVENT_TYPES = total_size(starting_states);
			std::vector<Resolvent<Matrix, Vector>> resolvents(N_RESOLVENT_TYPES);
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

		/**
		 * @brief Compute resolvents and collect residual eigenvector information.
		 *
		 * Similar to `compute_collective_modes()` but, instead of returning only the
		 * resolvent summaries, this version also computes Lanczos residuals and returns
		 * detailed residual eigenvector data for post-analysis.
		 *
		 * @tparam CheckHermitian If >0, enables runtime Hermiticity checks against precision 10^-CheckHermitian.
		 * @param n_lanczos_iterations Number of Lanczos iterations used for each resolvent.
		 * @return A pair where the first element is a vector of `ResolventReturnData`
		 *         and the second element is a `std::list` of `ResidualData` entries
		 *         containing eigenvectors and eigenvalues extracted from Lanczos residuals.
		 */
		template<int CheckHermitian = -1>
		std::pair<std::vector<ResolventReturnData>, std::list<ResidualData>> compute_collective_modes_with_residuals(unsigned int n_lanczos_iterations)
		{
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix, transform_matrix;

			std::pair<std::vector<ResolventReturnData>, std::list<ResidualData>> return_data;
			std::list<ResidualData>& residual_infos = return_data.second;

			set_begin();
			const int N_RESOLVENT_TYPES = total_size(starting_states);
			std::vector<Resolvent<Matrix, Vector>> resolvents(N_RESOLVENT_TYPES);
			auto resolvent_it = resolvents.begin();

			// No need to do anything if there are no starting states for the phase channel
			if (phase_size(starting_states) > 0U) {
				compute_solver_matrix<0, 1>(k_solutions, solver_matrix, transform_matrix);
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

			// No need to do anything if there are no starting states for the amplitude channel
			if (amplitude_size(starting_states) > 0U) {
				compute_solver_matrix<1, 0>(k_solutions, solver_matrix, transform_matrix);
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

		/**
		* @brief Performs a full diagonalization of the dynamical matrices.
		* Saves the first <n_residuals> eigenvectors, all eigenvalues and the weights for each starting state.
		* NOTE: This algorithm keeps only one eigenvalue (and its corresponding weight) from the nullspace of the solver matrix.
		* The nullspace is determined via the precision given to the constructor, see get_number_of_zero_eigenvalues().
		* @return A pair of FullDiagData objects, the first for the phase states, the second for the amplitude states.
		*/
		template<int CheckHermitian = -1>
		std::pair<FullDiagData, FullDiagData> full_diagonalization() {
			auto k_solutions = this->diagonalize_K_matrices<CheckHermitian>();
			Matrix solver_matrix, transform_matrix;

			Eigen::SelfAdjointEigenSolver<Matrix> solver;
			std::pair<FullDiagData, FullDiagData> return_data;

			// No need to do anything if there are no starting states for the phase channel
			if (phase_size(starting_states) > 0U) {
				compute_solver_matrix<0, 1>(k_solutions, solver_matrix, transform_matrix);
				set_begin();
				solver.compute(solver_matrix);
				print_duration("Time for first ED: ");
				TransformQR qr(transform_matrix); 
				print_duration("Time for first QR decomp: ");
				const std::size_t n_zero = get_number_of_zero_eigenvalues(solver);
				const std::size_t n_non_zero = solver.eigenvalues().size() - n_zero;
				return_data.first = set_full_diag_data<const_phase_it>(solver, qr, n_zero, n_non_zero, transform_matrix);
				if constexpr (check_qr) {
					for (std::size_t i = 0U; i < n_residuals; ++i) {
						Eigen::Map<Vector> _eigen(return_data.first.first_eigenvectors[i].data(), return_data.first.first_eigenvectors[i].size());
						const auto reconstructed = transform_matrix * _eigen;
						std::cout << "|solve_matrix * reconstructed eigenvector - eigenvalue * reconstructed eigenvector| = " 
							<< (solver_matrix * reconstructed - return_data.first.eigenvalues[i] * return_data.first.eigenvalues[i] * reconstructed).norm() << std::endl;
					}
				}
			
			}

			// No need to do anything if there are no starting states for the amplitude channel
			if (amplitude_size(starting_states) > 0U) {
				compute_solver_matrix<1, 0>(k_solutions, solver_matrix, transform_matrix);
				set_begin();
				solver.compute(solver_matrix);
				print_duration("Time for second ED: ");
				TransformQR qr(transform_matrix); 
				const std::size_t n_zero = get_number_of_zero_eigenvalues(solver);
				const std::size_t n_non_zero = solver.eigenvalues().size() - n_zero;

				print_duration("Time for second QR decomp: ");
				return_data.second = set_full_diag_data<const_amplitude_it>(solver, qr, n_zero, n_non_zero, transform_matrix);
				if constexpr (check_qr) {
					for (std::size_t i = 0U; i < n_residuals; ++i) {
						Eigen::Map<Vector> _eigen(return_data.second.first_eigenvectors[i].data(), return_data.second.first_eigenvectors[i].size());
						const auto reconstructed = transform_matrix * _eigen;
						std::cout << "|solve_matrix * reconstructed eigenvector - eigenvalue * reconstructed eigenvector| = " 
							<< (solver_matrix * reconstructed - return_data.second.eigenvalues[i] * return_data.second.eigenvalues[i] * reconstructed).norm() << std::endl;
					}
				}
			}

			return return_data;
		}
	};
}
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_XPRESOLVENT_HPP

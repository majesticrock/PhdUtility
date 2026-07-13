#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_RESOLVENT_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_RESOLVENT_HPP
// Use (void) to silence unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#include <type_traits>
#include <Eigen/Dense>
#include <cmath>


#include "detail/UnderlyingRealType.hpp"
#include "detail/is_complex.hpp"
#include "detail/GramSchmidt.hpp"
#include "ResolventDataTypes.hpp"

namespace mrock::iEoM {
	using std::abs;

	/**
	 * @brief Performs Lanczos-based resolvent computations for linear operators.
	 *
	 * This class builds recurrence coefficients for continued fraction resolvent
	 * representations. It supports Hermitian and symplectic problems and can
	 * store residual information for low-lying eigenvalues when requested.
	 *
	 * @tparam EigenMatrixType Matrix type used for the operator; must be square.
	 * @tparam EigenVectorType Vector type used for the starting state and Krylov basis.
	 */
	template <class EigenMatrixType, class EigenVectorType>
	class Resolvent
	{
	private:
		using RealType = detail::UnderlyingRealType_t<typename EigenMatrixType::Scalar>;
		using ComputationType = typename EigenMatrixType::Scalar;
		using resolvent_data = ResolventData<RealType>;
		static constexpr bool isComplex = detail::is_complex_v<typename EigenVectorType::Scalar>;

	public:
		/** Starting state vector used to generate the Lanczos basis. */
		EigenVectorType startingState;
		/** Computed resolvent data including recurrence coefficients and residual information. */
		ResolventDataWrapper<RealType> data;
		/**
		 * @brief Set the starting state used for the resolvent iteration.
		 *
		 * @param state New starting state vector.
		 */
		inline void set_starting_state(const EigenVectorType& state) {
			this->startingState = state;
		};
		/**
		 * @brief Get the current starting state.
		 *
		 * @return Const reference to the starting state vector.
		 */
		const EigenVectorType& getStartingState() const {
			return this->startingState;
		}
		/**
		 * @brief Construct with a starting state and optional name.
		 *
		 * @param _StargingState Initial Krylov vector for the resolvent iteration.
		 * @param name Optional name stored in the result wrapper.
		 */
		Resolvent(const EigenVectorType& _StargingState, const std::string& name="") 
			: startingState(_StargingState), data(name) { };
		/**
		 * @brief Construct with an optional name only.
		 *
		 * @param name Name stored in the result wrapper.
		 */
		Resolvent(const std::string name)
			: data(name) { };
		Resolvent() = default;

		/**
		 * @brief Compute generalized Lanczos coefficients for a symplectic problem.
		 *
		 * @param toSolve Operator matrix to solve.
		 * @param symplectic Symplectic metric matrix; must be positive semidefinite.
		 * @param maxIter Maximum Lanczos iterations.
		 */
		void compute(const EigenMatrixType& toSolve, const EigenMatrixType& symplectic, int maxIter)
		{
			const size_t matrix_size = toSolve.rows();
			maxIter = std::min(maxIter, static_cast<int>(matrix_size));

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			EigenVectorType currentSolution(matrix_size); // corresponds to |q_(i+1)>
			// First filling
			std::vector<EigenVectorType> basis_vectors;
			EigenVectorType first = EigenVectorType::Zero(matrix_size); // corresponds to |q_0>
			EigenVectorType second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			ComputationType norm_buffer = second.dot(symplectic * second);
			if constexpr (isComplex) {
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			}
			res.b_i.push_back(abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basis_vectors.push_back(first);
			basis_vectors.push_back(second);

			std::vector<RealType> alphas, betas;
			alphas.reserve(maxIter);
			betas.reserve(maxIter);

			betas.push_back(1);
			size_t iterNum{};
			bool goOn = true;
			EigenVectorType buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basis_vectors.back();
				norm_buffer = basis_vectors.back().dot(symplectic * buffer);
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
					alphas.push_back(norm_buffer.real());
				}
				else {
					alphas.push_back(norm_buffer);
				}

				currentSolution = (buffer - (alphas.back() * basis_vectors.back())) - (betas.back() * basis_vectors[iterNum]);
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				}
				betas.push_back(abs(norm_buffer));
				basis_vectors.push_back(currentSolution / betas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter || abs(betas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (long i = 0; i < alphas.size(); i++)
			{
				res.a_i.push_back(alphas[i]);
				res.b_i.push_back(betas[i + 1]);
			}
			data.push_back(std::move(res));
		}

		/**
		 * @brief Compute Lanczos coefficients for a Hermitian problem.
		 *
		 * This overload assumes the symplectic metric is the identity matrix.
		 *
		 * @param toSolve Operator matrix to solve.
		 * @param maxIter Maximum Lanczos iterations.
		 */
		void compute(const EigenMatrixType& toSolve, int maxIter)
		{
			const size_t matrix_size = toSolve.rows();
			maxIter = std::min(maxIter, static_cast<int>(matrix_size));

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			EigenVectorType currentSolution(matrix_size); // corresponds to |q_(i+1)>
			// First filling
			std::vector<EigenVectorType> basis_vectors;
			EigenVectorType first = EigenVectorType::Zero(matrix_size); // corresponds to |q_0>
			EigenVectorType second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

			second /= sqrt(res.b_i.back());
			basis_vectors.push_back(first);
			basis_vectors.push_back(second);

			std::vector<RealType> alphas, betas;
			alphas.reserve(maxIter);
			betas.reserve(maxIter);

			betas.push_back(1);
			size_t iterNum{};
			bool goOn = true;
			EigenVectorType buffer;

			while (goOn) {
				// algorithm
				buffer = toSolve * basis_vectors.back();
				if constexpr (isComplex) {
					// This has to be real, as <x|H|x> is always real if H=H^+
					alphas.push_back(basis_vectors.back().dot(buffer).real());
				}
				else {
					alphas.push_back(basis_vectors.back().dot(buffer));
				}
				currentSolution = buffer - (alphas.back() * basis_vectors.back() + betas.back() * basis_vectors[iterNum]);
				betas.push_back(currentSolution.norm());
				basis_vectors.push_back(currentSolution / betas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter || abs(betas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (size_t i = 0U; i < alphas.size(); ++i)
			{
				res.a_i.push_back(alphas[i]);
				res.b_i.push_back(betas[i + 1] * betas[i + 1]);
			}
			// The last b is irrelevant, it does not really exist; it's an artifact of the algorithm
			res.b_i.pop_back();
			data.push_back(std::move(res));
		}

		/**
		 * @brief Compute generalized Lanczos coefficients using an explicit M and N representation.
		 *
		 * This form can be more stable for complex-valued matrices whose inner product is defined by N.
		 *
		 * @param toSolve Operator matrix to solve.
		 * @param symplectic Symplectic metric matrix; must be positive semidefinite.
		 * @param N Auxiliary matrix used in the formulation.
		 * @param maxIter Maximum Lanczos iterations.
		 */
		void compute_from_N_M(const EigenMatrixType& toSolve, const EigenMatrixType& symplectic, const EigenMatrixType& N, int maxIter)
		{
			auto matrix_size = toSolve.rows();
			maxIter = std::min(maxIter, static_cast<int>(matrix_size));

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			EigenVectorType currentSolution(matrix_size); // corresponds to |q_(i+1)>
			// First filling
			std::vector<EigenVectorType> basis_vectors;
			EigenVectorType first = EigenVectorType::Zero(matrix_size); // corresponds to |q_0>
			EigenVectorType second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			ComputationType norm_buffer = second.dot(symplectic * second);
			if constexpr (isComplex) {
				assertm(abs(norm_buffer.imag()) < 1e-6, "First norm is complex! ");
			}
			res.b_i.push_back(abs(norm_buffer));

			second /= sqrt(res.b_i.back());
			basis_vectors.push_back(first);
			basis_vectors.push_back(second);

			std::vector<RealType> alphas, betas;
			alphas.reserve(maxIter);
			betas.reserve(maxIter);

			betas.push_back(1);
			size_t  iterNum{};
			bool goOn = true;
			EigenVectorType buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basis_vectors.back();
				norm_buffer = basis_vectors.back().dot(N * basis_vectors.back());
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "First norm in loop is complex!");
					alphas.push_back(norm_buffer.real());
				}
				else {
					alphas.push_back(norm_buffer);
				}

				currentSolution = (buffer - (alphas.back() * basis_vectors.back())) - (betas.back() * basis_vectors[iterNum]);
				norm_buffer = sqrt(currentSolution.dot(symplectic * currentSolution));
				if constexpr (isComplex) {
					assertm(abs(norm_buffer.imag()) < 1e-6, "Second norm in loop is complex!");
				}
				betas.push_back(abs(norm_buffer));
				basis_vectors.push_back(currentSolution / betas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter || abs(betas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (size_t i = 0U; i < alphas.size(); ++i)
			{
				res.a_i.push_back(alphas[i]);
				res.b_i.push_back(betas[i + 1]);
			}
			data.push_back(std::move(res));
		}

		/**
		 * @brief Compute Lanczos coefficients for a Hermitian problem with basis reorthogonalization.
		 *
		 * This version reorthogonalizes the Krylov basis at every iteration to improve numerical stability.
		 *
		 * @param toSolve Hermitian operator matrix.
		 * @param maxIter Maximum Lanczos iterations.
		 */
		void compute_with_reorthogonalization(const EigenMatrixType& toSolve, int maxIter)
		{
			const size_t matrix_size = toSolve.rows();
			maxIter = std::min(maxIter, static_cast<int>(matrix_size));

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			EigenVectorType currentSolution(matrix_size); // corresponds to |q_(i+1)>
			// First filling
			std::vector<EigenVectorType> basis_vectors;
			EigenVectorType second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

			second /= sqrt(res.b_i.back());
			basis_vectors.push_back(second);

			std::vector<RealType> alphas, betas;
			alphas.reserve(maxIter);
			betas.reserve(maxIter);

			betas.push_back(1);
			size_t iterNum{};
			bool goOn = true;
			EigenVectorType buffer;
			while (goOn) {
				// algorithm
				buffer = toSolve * basis_vectors.back();
				if constexpr (isComplex) {
					// This has to be real, as <x|H|x> is always real if H=H^+
					alphas.push_back(basis_vectors.back().dot(buffer).real());
				}
				else {
					alphas.push_back(basis_vectors.back().dot(buffer));
				}
				if (iterNum > 0U) {
					currentSolution = buffer - (alphas.back() * basis_vectors.back() + betas.back() * basis_vectors[iterNum]);
				}
				else {
					currentSolution = buffer - (alphas.back() * basis_vectors.back());
				}
				detail::GramSchmidt<typename EigenVectorType::Scalar>::orthogonalize_single_vector(currentSolution, basis_vectors);
				betas.push_back(currentSolution.norm());
				basis_vectors.push_back(currentSolution / betas.back());
				++iterNum;

				// breaking conditions
				if (iterNum >= maxIter || abs(betas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (size_t i = 0U; i < alphas.size(); ++i)
			{
				res.a_i.push_back(alphas[i]);
				res.b_i.push_back(betas[i + 1] * betas[i + 1]);
			}
			// The last b is irrelevant, it does not really exist; it's an artifact of the algorithm
			res.b_i.pop_back();
			data.push_back(std::move(res));
		}

		/**
		 * @brief Compute Lanczos coefficients and low-lying residuals for a Hermitian problem.
		 *
		 * This method uses reorthogonalization and periodically estimates residuals for the
		 * lowest @p n_residuals Ritz values.
		 *
		 * @tparam n_residuals Number of residual eigenvalues to track.
		 * @param toSolve Hermitian operator matrix.
		 * @param maxIter Maximum Lanczos iterations.
		 * @return ResidualInformation containing eigenvalue residuals, converged eigenvectors, and weights.
		 */
		template <int n_residuals>
		ResidualInformation<RealType, n_residuals> compute_with_residuals(const EigenMatrixType& toSolve, int maxIter)
		{
			static_assert(isComplex == false, "Residual computation not implemented for complex types yet!");
			const size_t matrix_size = toSolve.rows();
			maxIter = std::min(maxIter, static_cast<int>(matrix_size));

			if (toSolve.rows() != toSolve.cols()) {
				std::cerr << "Matrix is not square!" << std::endl;
				throw;
			}

			EigenVectorType currentSolution(matrix_size); // corresponds to |q_(i+1)>
			// First filling
			std::vector<EigenVectorType> basis_vectors;
			EigenVectorType second = this->startingState; // corresponds to |q_1>
			resolvent_data res;
			res.b_i.push_back(second.squaredNorm());

			second /= sqrt(res.b_i.back());
			basis_vectors.push_back(second);

			std::vector<RealType> alphas, betas;
			alphas.reserve(maxIter);
			betas.reserve(maxIter);

			betas.push_back(1);
			int iterNum{};
			bool goOn = true;
			EigenVectorType buffer;

			Eigen::SelfAdjointEigenSolver<EigenMatrixType> eigen_solver;
			typedef Eigen::Vector<RealType, Eigen::Dynamic> TVector;
			typedef Eigen::Map<const TVector> TMap;

			ResidualInformation<RealType, n_residuals> residual_info;

			while (goOn) {
				// algorithm
				buffer = toSolve * basis_vectors.back();
				if constexpr (isComplex) {
					// This has to be real, as <x|H|x> is always real if H=H^+
					alphas.push_back(basis_vectors.back().dot(buffer).real());
				}
				else {
					alphas.push_back(basis_vectors.back().dot(buffer));
				}
				if (iterNum > 0U) {
					currentSolution = buffer - (alphas.back() * basis_vectors.back() + betas.back() * basis_vectors[iterNum]);
				}
				else {
					currentSolution = buffer - (alphas.back() * basis_vectors.back());
				}

				detail::GramSchmidt<typename EigenVectorType::Scalar>::orthogonalize_single_vector(currentSolution, basis_vectors);
				betas.push_back(currentSolution.norm());
				basis_vectors.push_back(currentSolution / betas.back());
				++iterNum;

				if (iterNum > n_residuals + 1 && iterNum % 10 == 0) {
					if (!std::all_of(residual_info.converged.begin(), residual_info.converged.end(), [](bool v) { return v; })) {
						eigen_solver.computeFromTridiagonal(TMap(alphas.data(), iterNum), TMap(betas.data() + 1, iterNum - 1));
						int skip{};
						int i_skip{};
						// see https://epubs.siam.org/doi/book/10.1137/1.9780898719581, chapter 4.4, eq. 4.13
						for (int i = 0; i < n_residuals; ++i) {
							if (residual_info.converged[i]) continue;

							redo_if_ghost:
							i_skip = i + skip; // To avoid Lanczos ghosts
							if (i_skip >= iterNum) break; 
							residual_info.residuals[i] = betas.back() * abs(eigen_solver.eigenvectors().col(i_skip)(iterNum - 1));
							constexpr double CONVERGE_EPS = 1.49011612e-8; // ~ sqrt(machine epsilon) ~ 2^{-26}
							if (residual_info.residuals[i] < CONVERGE_EPS || iterNum >= maxIter) {
								// Eigen sorts the eigenvalues in ascending order
								for (int c = 0; c < n_residuals; ++c) { //&& residual_info.converged[c]
									if (abs(residual_info.eigenvalues[c] - eigen_solver.eigenvalues()(i_skip)) < 1e-12) {
										// Found a Lanczos ghost
										++skip;
										goto redo_if_ghost;
									}
								}

								residual_info.converged[i] = residual_info.residuals[i] < CONVERGE_EPS;
								residual_info.eigenvalues[i] = eigen_solver.eigenvalues()(i_skip);
								residual_info.n_ghosts[i] = skip;

								// Computing the eigenvector in the original space
								EigenVectorType eigvec = EigenVectorType::Zero(matrix_size);
								for (int j = 0; j < iterNum; ++j) {
									eigvec += eigen_solver.eigenvectors().col(i_skip)(j) * basis_vectors[j];
								}
								// The eigenvector is useless, if the residual is not small
								// However, the eigenvalue is bounded by lambda_true = lambda_approx +/- residual
								if (residual_info.converged[i]) { 
									residual_info.eigenvectors[i] = std::vector<RealType>(eigvec.data(), eigvec.data() + eigvec.size());
									residual_info.weights[i] = eigvec.dot(this->startingState);
								}
							}
							if (!residual_info.converged[i] && iterNum < maxIter) break; // Fill list form the bottom up
						}
					}
				}

				// breaking conditions
				if (iterNum >= maxIter || abs(betas.back()) < 1e-10) {
					goOn = false;
				}
			}
			for (size_t i = 0U; i < alphas.size(); ++i)
			{
				res.a_i.push_back(alphas[i]);
				res.b_i.push_back(betas[i + 1] * betas[i + 1]);
			}
			// The last b is irrelevant, it does not really exist; it's an artifact of the algorithm
			res.b_i.pop_back();
			data.push_back(std::move(res));

			return residual_info;
		}

		/**
		 * @brief Get the computed resolvent data.
		 *
		 * @return Const reference to the stored ResolventDataWrapper containing
		 *         continued fraction coefficients and residual information.
		 */
		const ResolventDataWrapper<RealType>& get_data() const {
			return data;
		}
	};
}
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_RESOLVENT_HPP

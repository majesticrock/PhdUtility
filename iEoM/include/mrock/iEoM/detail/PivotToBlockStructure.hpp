#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_PIVOTTOBLOCKSTRUCTURE_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_PIVOTTOBLOCKSTRUCTURE_HPP
#include <Eigen/Dense>

#include "UnderlyingRealType.hpp"
#include "BlockDiagonalMatrix.hpp"

namespace mrock::iEoM::detail {
	/**
	 * @brief Compute a permutation that groups zero off-diagonal blocks.
	 *
	 * The returned permutation reorders rows and columns such that
	 * contiguous zero off-diagonal regions in the matrix become block
	 * diagonal structure.
	 *
	 * @tparam EigenMatrixType Matrix type to analyze.
	 * @param matrix Input matrix.
	 * @param epsilon Tolerance below which entries are treated as zero.
	 * @return Permutation matrix representing the pivot.
	 */
	template<class EigenMatrixType>
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> pivot_to_block_structure(const EigenMatrixType& matrix, const detail::RealScalar<EigenMatrixType> epsilon = 1e-12) {
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P(matrix.rows());
		P.setIdentity();

		auto permuted_matrix_element = [&matrix, &P](Eigen::Index i, Eigen::Index j) {
			return matrix.coeffRef(P.indices()(i), P.indices()(j));
			};

		for (int i = 0; i < matrix.rows(); ++i) {
			int offset = 1;
			for (int j = i + 1; j < matrix.cols(); ++j) {
				if (abs(permuted_matrix_element(i, j)) > epsilon) {
					P.applyTranspositionOnTheRight(j, i + offset);
					++offset;
				}
			}
		}
		return P;
	};


	/**
	 * @brief Wrapper for eigenvector and eigenvalue results.
	 *
	 * Provides a common interface for solving both dense and block
	 * diagonal Hermitian matrices.
	 *
	 * @tparam MatrixType Matrix type used for eigenvectors.
	 * @tparam RealType Underlying real scalar type for eigenvalues.
	 * @tparam RealVector Vector type used to store eigenvalues.
	 */
	template <class MatrixType, class RealType = UnderlyingRealType_t<typename MatrixType::Scalar>, class RealVector = Eigen::Vector<RealType, Eigen::Dynamic>>
	struct matrix_wrapper {
		MatrixType eigenvectors;
		RealVector eigenvalues;

		/**
		 * @brief Default-construct an empty result wrapper.
		 */
		inline matrix_wrapper() {};

		/**
		 * @brief Construct a result wrapper with preallocated storage.
		 *
		 * @param size Dimension of the matrix to reconstruct or solve.
		 */
		inline explicit matrix_wrapper(Eigen::Index size)
			: eigenvectors(MatrixType::Zero(size, size)), eigenvalues(RealVector::Zero(size))
		{};

		/**
		 * @brief Reconstruct the full matrix from eigenvectors and eigenvalues.
		 *
		 * @return Reconstructed matrix equal to V * D * V^H.
		 */
		inline MatrixType reconstruct_matrix() const
		{
			return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.adjoint();
		};

		/**
		 * @brief Solve each Hermitian block independently.
		 *
		 * @param toSolve Full matrix to solve, containing block diagonal structure.
		 * @param blocks Identified Hermitian blocks within the matrix.
		 * @return matrix_wrapper containing eigenpairs for the full matrix.
		 */
		static matrix_wrapper solve_block_diagonal_matrix(const MatrixType& toSolve, const std::vector<HermitianBlock>& blocks) {
			matrix_wrapper solution(toSolve.rows());

#ifdef MROCK_IEOM_PARALLELIZE_BLOCKMATRIX
#pragma omp parallel for
#endif
			for (int i = 0; i < blocks.size(); ++i)
			{
				Eigen::SelfAdjointEigenSolver<MatrixType> solver(toSolve.block(blocks[i].position, blocks[i].position, blocks[i].size, blocks[i].size));
				solution.eigenvalues.segment(blocks[i].position, blocks[i].size) = solver.eigenvalues();
				solution.eigenvectors.block(blocks[i].position, blocks[i].position, blocks[i].size, blocks[i].size) = solver.eigenvectors();
			}
			
			return solution;
		}

		/**
		 * @brief Pivot the matrix, solve each block, and restore eigenvectors.
		 *
		 * The input matrix is permuted to block diagonal form, solved block-wise,
		 * and the eigenvector basis is transformed back to the original ordering.
		 *
		 * @param toSolve Matrix to solve; modified in-place.
		 * @return matrix_wrapper containing the eigenpairs of the original matrix.
		 */
		static matrix_wrapper pivot_and_solve(MatrixType& toSolve)
		{
			auto pivot = pivot_to_block_structure(toSolve);
			toSolve = pivot.transpose() * toSolve * pivot;
			const auto blocks = identify_hermitian_blocks(toSolve);
			auto solution = solve_block_diagonal_matrix(toSolve, blocks);
			solution.eigenvectors.applyOnTheLeft(pivot);
			return solution;
		};

		/**
		 * @brief Solve the full Hermitian matrix without pivoting.
		 *
		 * @param toSolve Matrix to diagonalize.
		 * @return matrix_wrapper containing eigenpairs.
		 */
		static matrix_wrapper only_solve(MatrixType& toSolve)
		{
			Eigen::SelfAdjointEigenSolver<MatrixType> solver(toSolve);
			matrix_wrapper solution(toSolve.rows());
			solution.eigenvalues = solver.eigenvalues();
			solution.eigenvectors = solver.eigenvectors();
			return solution;
		};

		/**
		 * @brief Test whether the matrix is non-negative definite.
		 *
		 * The matrix is pivoted to block diagonal form, and each block's
		 * smallest eigenvalue is compared against -EPSILON.
		 *
		 * @param toSolve Matrix to test.
		 * @param EPSILON Tolerance for negative eigenvalues.
		 * @return True if the matrix is non-negative definite.
		 */
		static bool is_non_negative(MatrixType& toSolve, const RealType EPSILON)
		{
			auto pivot = pivot_to_block_structure(toSolve);
			toSolve = pivot.transpose() * toSolve * pivot;
			auto blocks = identify_hermitian_blocks(toSolve);
			for (const auto& block : blocks)
			{
				Eigen::SelfAdjointEigenSolver<MatrixType> solver(toSolve.block(block.position, block.position, block.size, block.size), Eigen::EigenvaluesOnly);

				if ((solver.eigenvalues().array() < -EPSILON).any()) {
					return false;
				}
			}
			return true;
		};
	};

	/**
	 * @brief Specialization of matrix_wrapper for block diagonal matrices.
	 *
	 * This wrapper stores eigenvectors using a BlockDiagonalMatrix representation
	 * and supports reconstruction back to a dense matrix if required.
	 *
	 * @tparam Number Scalar type of the underlying matrix entries.
	 */
	template <class Number>
	struct matrix_wrapper<BlockDiagonalMatrix<detail::MatrixN<Number>>, UnderlyingRealType_t<Number>, Eigen::Vector<UnderlyingRealType_t<Number>, Eigen::Dynamic>> {
		BlockDiagonalMatrix<detail::MatrixN<Number>> eigenvectors;
		Eigen::Vector<UnderlyingRealType_t<Number>, Eigen::Dynamic> eigenvalues;

		/**
		 * @brief Reconstruct the matrix as a block diagonal object.
		 *
		 * @return BlockDiagonalMatrix representing the reconstructed matrix.
		 */
		inline auto reconstruct_matrix() const
		{
			return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.adjoint();
		}
		/**
		 * @brief Reconstruct the full dense matrix from the block diagonal form.
		 *
		 * @return Dense matrix assembled from the block diagonal eigenstructure.
		 */
		inline detail::MatrixN<Number> reconstruct_matrix_as_eigen() const {
			return reconstruct_matrix().construct_matrix();
		}

		/**
		 * @brief Solve a block diagonal matrix by computing eigenpairs block-wise.
		 *
		 * @param toSolve Block diagonal matrix to solve.
		 * @return matrix_wrapper containing block eigenvectors and eigenvalues.
		 */
		static matrix_wrapper solve_block_diagonal_matrix(const BlockDiagonalMatrix<detail::MatrixN<Number>>& toSolve) {
			matrix_wrapper solution;
			solution.eigenvalues = Eigen::Vector<UnderlyingRealType_t<Number>, Eigen::Dynamic>::Zero(toSolve.rows());
			solution.eigenvectors.blocks.resize(toSolve.blocks.size());
			solution.eigenvectors.blocks_begin = toSolve.blocks_begin;
#ifdef MROCK_IEOM_PARALLELIZE_BLOCKMATRIX
#pragma omp parallel for
#endif
			for (int i = 0; i < toSolve.blocks.size(); ++i)
			{
				Eigen::SelfAdjointEigenSolver<detail::MatrixN<Number>> solver(toSolve.blocks[i]);
				solution.eigenvectors.blocks[i] = solver.eigenvectors();
				solution.eigenvalues.segment(toSolve.blocks_begin[i], toSolve.blocks[i].rows()) = solver.eigenvalues();
			}
			
			return solution;
		}
	};

	template <class Number> 
	using blocked_matrix_wrapper = matrix_wrapper<BlockDiagonalMatrix<detail::MatrixN<Number>>, UnderlyingRealType_t<Number>, Eigen::Vector<UnderlyingRealType_t<Number>, Eigen::Dynamic>>;
}
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_PIVOTTOBLOCKSTRUCTURE_HPP

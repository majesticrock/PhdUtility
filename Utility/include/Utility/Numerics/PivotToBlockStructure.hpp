#pragma once
#include <Eigen/Dense>
#include <omp.h>
#include "BlockDiagonalMatrix.hpp"

namespace Utility::Numerics {
	// Pivots a matrix so that all offdiagonal 0 blocks are contiguous
	// The permutation matrix is returned
	// epsilon is used to determine if a matrix element is 0 (especially important for floating point operations)
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


	template <class MatrixType, class RealType = UnderlyingFloatingPoint_t<typename MatrixType::Scalar>, class RealVector = Eigen::Vector<RealType, Eigen::Dynamic>>
	struct matrix_wrapper {
		MatrixType eigenvectors;
		RealVector eigenvalues;

		inline matrix_wrapper() {};

		inline explicit matrix_wrapper(Eigen::Index size)
			: eigenvectors(MatrixType::Zero(size, size)), eigenvalues(RealVector::Zero(size))
		{};

		inline MatrixType reconstruct_matrix() const
		{
			return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.adjoint();
		};

		static matrix_wrapper solve_block_diagonal_matrix(const MatrixType& toSolve, const std::vector<HermitianBlock>& blocks) {
			matrix_wrapper solution(toSolve.rows());

#ifdef _BLOCKS_USE_OMP
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

		static matrix_wrapper pivot_and_solve(MatrixType& toSolve)
		{
			auto pivot = pivot_to_block_structure(toSolve);
			toSolve = pivot.transpose() * toSolve * pivot;
			const auto blocks = identify_hermitian_blocks(toSolve);
			auto solution = solve_block_diagonal_matrix(toSolve, blocks);
			solution.eigenvectors.applyOnTheLeft(pivot);
			return solution;
		};

		static matrix_wrapper only_solve(MatrixType& toSolve)
		{
			Eigen::SelfAdjointEigenSolver<MatrixType> solver(toSolve);
			matrix_wrapper solution(toSolve.rows());
			solution.eigenvalues = solver.eigenvalues();
			solution.eigenvectors = solver.eigenvectors();
			return solution;
		};

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

	template <class Number>
	struct matrix_wrapper<BlockDiagonalMatrix<Number>, UnderlyingFloatingPoint_t<Number>, Eigen::Vector<UnderlyingFloatingPoint_t<Number>, Eigen::Dynamic>> {
		BlockDiagonalMatrix<Number> eigenvectors;
		Eigen::Vector<UnderlyingFloatingPoint_t<Number>, Eigen::Dynamic> eigenvalues;

		inline BlockDiagonalMatrix<Number> reconstruct_matrix() const
		{
			return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.adjoint();
		}
		inline BlockDiagonalMatrix<Number>::InternalMatrix reconstruct_matrix_as_eigen() const {
			return reconstruct_matrix().construct_matrix();
		}

		static matrix_wrapper solve_block_diagonal_matrix(const BlockDiagonalMatrix<Number>& toSolve) {
			matrix_wrapper solution;
			solution.eigenvalues = Eigen::Vector<UnderlyingFloatingPoint_t<Number>, Eigen::Dynamic>::Zero(toSolve.rows());
			solution.eigenvectors.blocks.resize(toSolve.blocks.size());
			solution.eigenvectors.blocks_begin = toSolve.blocks_begin;
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < toSolve.blocks.size(); ++i)
			{
				Eigen::SelfAdjointEigenSolver<typename BlockDiagonalMatrix<Number>::InternalMatrix> solver(toSolve.blocks[i]);
				solution.eigenvectors.blocks[i] = solver.eigenvectors();
				solution.eigenvalues.segment(toSolve.blocks_begin[i], toSolve.blocks[i].rows()) = solver.eigenvalues();
			}
			
			return solution;
		}
	};

	template <class Number> using blocked_matrix_wrapper = matrix_wrapper<BlockDiagonalMatrix<Number>, UnderlyingFloatingPoint_t<Number>, Eigen::Vector<UnderlyingFloatingPoint_t<Number>, Eigen::Dynamic>>;
}
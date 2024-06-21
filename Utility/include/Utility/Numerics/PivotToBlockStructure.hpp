#pragma once
#include <Eigen/Dense>
#include <omp.h>
#include "../UnderlyingFloatingPoint.hpp"

namespace Utility::Numerics {
	namespace Detail {
		template<class EigenMatrixType>
		using RealScalar = UnderlyingFloatingPoint_t<typename EigenMatrixType::Scalar>;
	}

	// Pivots a matrix so that all offdiagonal 0 blocks are contiguous
	// The permutation matrix is returned
	// epsilon is used to determine if a matrix element is 0 (especially important for floating point operations)
	template<class EigenMatrixType>
	auto pivot_to_block_structure(const EigenMatrixType& matrix, const Detail::RealScalar<EigenMatrixType> epsilon = 1e-12) {
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

	struct HermitianBlock {
		Eigen::Index position{};
		Eigen::Index size{};
	};

	inline std::ostream& operator<<(std::ostream& os, const HermitianBlock& block) {
		os << block.position << "\t" << block.size;
		return os;
	};

	template<class EigenMatrixType>
	auto identify_hermitian_blocks(const EigenMatrixType& matrix, const Detail::RealScalar<EigenMatrixType> epsilon = 1e-12) {
		Eigen::Index block_index{ 1 };
		Eigen::Index block_size{};
		std::vector<HermitianBlock> block_indices;
		for (Eigen::Index i = 0; i < matrix.rows(); ++i)
		{
			for (Eigen::Index j = matrix.cols() - 1; j > block_index; --j)
			{
				if (abs(matrix(i, j)) > epsilon) {
					block_index = j;
					break;
				}
			}
			if (block_index == i) {
				if (block_indices.empty()) {
					block_indices.push_back({ Eigen::Index{}, block_index + 1 });
				}
				else {
					block_size = block_index - (block_indices.back().size + block_indices.back().position) + 1;
					block_indices.push_back({ block_index - block_size + 1, block_size });
				}
			}
		}
		return block_indices;
	};

	template <class NumberType>
	struct matrix_wrapper {
		using RealType = UnderlyingFloatingPoint_t<NumberType>;
		using MatrixType = Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>;
		using RealVector = Eigen::Vector<RealType, Eigen::Dynamic>;
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

		static matrix_wrapper pivot_and_solve(MatrixType& toSolve)
		{
			auto pivot = pivot_to_block_structure(toSolve);
			toSolve = pivot.transpose() * toSolve * pivot;
			auto blocks = identify_hermitian_blocks(toSolve);
			matrix_wrapper solution(toSolve.rows());

#pragma omp parallel for
			for (int i = 0; i < blocks.size(); ++i)
			{
				Eigen::SelfAdjointEigenSolver<MatrixType> solver(toSolve.block(blocks[i].position, blocks[i].position, blocks[i].size, blocks[i].size));
				solution.eigenvalues.segment(blocks[i].position, blocks[i].size) = solver.eigenvalues();
				solution.eigenvectors.block(blocks[i].position, blocks[i].position, blocks[i].size, blocks[i].size) = solver.eigenvectors();
			}
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
}
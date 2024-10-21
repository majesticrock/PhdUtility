#pragma once
#include <Eigen/Dense>
#include <vector>
#include "../UnderlyingFloatingPoint.hpp"

namespace Utility::Numerics {
    namespace detail {
        template<class Number> using MatrixN = Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>;
        template<class Number> using VectorN = Eigen::Vector<Number, Eigen::Dynamic>;

        template<class EigenMatrixType> using RealScalar = UnderlyingFloatingPoint_t<typename EigenMatrixType::Scalar>;
    }

    struct HermitianBlock {
		Eigen::Index position{};
		Eigen::Index size{};
	};

    inline std::ostream& operator<<(std::ostream& os, const HermitianBlock& block) {
		os << block.position << "\t" << block.size;
		return os;
	};

	template<class EigenMatrixType>
	std::vector<HermitianBlock> identify_hermitian_blocks(const EigenMatrixType& matrix, const detail::RealScalar<EigenMatrixType> epsilon = 1e-12) {
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

    template <class Number>
	struct BlockDiagonalMatrix {
        using InternalMatrix = detail::MatrixN<Number>;
        using Scalar = Number;

        std::vector<InternalMatrix> blocks;
        std::vector<Eigen::Index> blocks_begin;

        BlockDiagonalMatrix() = default;

        BlockDiagonalMatrix(const detail::MatrixN<Number>& matrix, const std::vector<HermitianBlock>& block_indizes) {
            blocks.reserve(block_indizes.size());
            blocks_begin.reserve(block_indizes.size());
            for (const auto& block_index : block_indizes)
            {
                blocks.push_back(matrix.block(block_index.position, block_index.position, block_index.size, block_index.size));
                blocks_begin.push_back(block_index.position);
            }
        }
        BlockDiagonalMatrix(const detail::MatrixN<Number>& matrix) {
            const std::vector<HermitianBlock>& block_indizes = identify_hermitian_blocks(matrix);
            blocks.reserve(block_indizes.size());
            blocks_begin.reserve(block_indizes.size());
            for (const auto& block_index : block_indizes)
            {
                blocks.push_back(matrix.block(block_index.position, block_index.position, block_index.size, block_index.size));
                blocks_begin.push_back(block_index.position);
            }
        }

        InternalMatrix construct_matrix() const {
            InternalMatrix ret = InternalMatrix::Zero(rows(), cols());
            for(size_t i = 0U; i < blocks.size(); ++i) {
                ret.block(blocks_begin[i], blocks_begin[i], blocks[i].rows(), blocks[i].cols()) = blocks[i];
            }
            return ret;
        }

        inline Eigen::Index rows() const {
            Eigen::Index __rows{};
            for(const auto& block : blocks) {
                __rows += block.rows();
            }
            return __rows;
        }
        inline Eigen::Index cols() const {
            Eigen::Index __cols{};
            for(const auto& block : blocks) {
                __cols += block.cols();
            }
            return __cols;
        }

        BlockDiagonalMatrix& operator+=(const BlockDiagonalMatrix& rhs) {
            for(int i = 0U; i < blocks.size(); ++i) {
                this->blocks[i] += rhs.blocks[i];
            }
            return *this;
        }
        BlockDiagonalMatrix& operator-=(const BlockDiagonalMatrix& rhs) {
            for(int i = 0U; i < blocks.size(); ++i) {
                this->blocks[i] -= rhs.blocks[i];
            }
            return *this;
        }
        BlockDiagonalMatrix& operator*=(const BlockDiagonalMatrix& rhs) {
            for(int i = 0U; i < blocks.size(); ++i) {
                this->blocks[i] *= rhs.blocks[i];
            }
            return *this;
        }
	};

    template <class Number>
    BlockDiagonalMatrix<Number> operator+(BlockDiagonalMatrix<Number> lhs, const BlockDiagonalMatrix<Number>& rhs) {
        return (lhs += rhs);
    }
    template <class Number>
    BlockDiagonalMatrix<Number> operator-(BlockDiagonalMatrix<Number> lhs, const BlockDiagonalMatrix<Number>& rhs) {
        return (lhs -= rhs);
    }
    template <class Number>
    BlockDiagonalMatrix<Number> operator*(BlockDiagonalMatrix<Number> lhs, const BlockDiagonalMatrix<Number>& rhs) {
        return (lhs *= rhs);
    }

    template <class Number>
    detail::MatrixN<Number> operator*(detail::MatrixN<Number> lhs, const BlockDiagonalMatrix<detail::MatrixN<Number>>& rhs) {
        for(size_t i = 0U; i < rhs.blocks.size(); ++i){
            lhs.block(rhs.blocks_begin[i], rhs.blocks_begin[i], lhs.rows(), rhs.blocks[i].size()).applyOnTheRight(rhs.blocks[i]);
        }
        return lhs;
    }

    template <class Number>
    detail::MatrixN<Number> operator*(const BlockDiagonalMatrix<detail::MatrixN<Number>>& lhs, detail::MatrixN<Number> rhs) {
        for(size_t i = 0U; i < lhs.blocks.size(); ++i){
            rhs.block(lhs.blocks_begin[i], lhs.blocks_begin[i], lhs.blocks[i].size(), rhs.cols()).applyOnTheLeft(lhs.blocks[i]);
        }
        return rhs;
    }
}
#pragma once
#include <Eigen/Dense>
#include <vector>
#include "../UnderlyingFloatingPoint.hpp"

namespace Utility::Numerics {
    namespace detail {
        template<class Number> using MatrixN = Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>;

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

    template <class Number> struct BlockDiagonalMatrix;
    template <class Number>
    struct BlockDiagonalMatrixAdjointView {
        using InternalMatrix = detail::MatrixN<Number>;
        using Scalar = Number;

        const BlockDiagonalMatrix<Number>& _matrix;

        BlockDiagonalMatrixAdjointView() = default;
        explicit BlockDiagonalMatrixAdjointView(const BlockDiagonalMatrix<Number>& matrix) : _matrix(matrix) {};

        InternalMatrix construct_matrix() const {
            return _matrix.construct_matrix.adjoint();
        }
        inline Eigen::Index rows() const {
            return _matrix.cols();
        }
        inline Eigen::Index cols() const {
            return _matrix.rows();
        }
	};

    template <class Number>
	struct BlockDiagonalMatrix {
        using InternalMatrix = detail::MatrixN<Number>;
        using Scalar = Number;

        std::vector<InternalMatrix> blocks;
        std::vector<Eigen::Index> blocks_begin;

        BlockDiagonalMatrix() = default;
        BlockDiagonalMatrix(const detail::MatrixN<Number>& matrix, const std::vector<HermitianBlock>& block_indizes) 
        {
            blocks.reserve(block_indizes.size());
            blocks_begin.reserve(block_indizes.size());
            for (const auto& block_index : block_indizes)
            {
                blocks.push_back(matrix.block(block_index.position, block_index.position, block_index.size, block_index.size));
                blocks_begin.push_back(block_index.position);
            }
        }
        explicit BlockDiagonalMatrix(const detail::MatrixN<Number>& matrix) 
        {
            const std::vector<HermitianBlock>& block_indizes = identify_hermitian_blocks(matrix);
            blocks.reserve(block_indizes.size());
            blocks_begin.reserve(block_indizes.size());
            for (const auto& block_index : block_indizes)
            {
                blocks.push_back(matrix.block(block_index.position, block_index.position, block_index.size, block_index.size));
                blocks_begin.push_back(block_index.position);
            }
        }
        BlockDiagonalMatrix(const BlockDiagonalMatrixAdjointView<Number>& _adjoint) 
            : blocks_begin(_adjoint._matrix.blocks_begin) 
        {
            this->blocks.reserve(_adjoint._matrix.blocks.size());
            for(const auto& block : _adjoint._matrix.blocks) {
                this->blocks.push_back(block.adjoint());
            }
        }

        InternalMatrix construct_matrix() const 
        {
            InternalMatrix ret = InternalMatrix::Zero(rows(), cols());
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                ret.block(blocks_begin[i], blocks_begin[i], blocks[i].rows(), blocks[i].cols()) = blocks[i];
            }
            return ret;
        }

        BlockDiagonalMatrixAdjointView<Number> adjoint() const {
            return BlockDiagonalMatrixAdjointView<Number>(*this);
        }
        BlockDiagonalMatrix& adjointInPlace() {
            for(auto& block : this->blocks) {
                block.adjointInPlace();
            }
            return *this;
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

        BlockDiagonalMatrix& applyOnTheLeft(BlockDiagonalMatrix other) {
            other *= (*this);
            return (*this = other);
        }

        BlockDiagonalMatrix operator-() const {
            BlockDiagonalMatrix copy(*this);
            for(auto& block : copy.blocks) {
                block *= -1;
            }
            return copy;
        }

        BlockDiagonalMatrix& operator+=(const BlockDiagonalMatrix& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] += rhs.blocks[i];
            }
            return *this;
        }
        BlockDiagonalMatrix& operator-=(const BlockDiagonalMatrix& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] -= rhs.blocks[i];
            }
            return *this;
        }
        BlockDiagonalMatrix& operator*=(const BlockDiagonalMatrix& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] *= rhs.blocks[i];
            }
            return *this;
        }
	
        BlockDiagonalMatrix& operator+=(const BlockDiagonalMatrixAdjointView<Number>& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] += rhs._matrix.blocks[i].adjoint();
            }
            return *this;
        }
        BlockDiagonalMatrix& operator-=(const BlockDiagonalMatrixAdjointView<Number>& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] -= rhs._matrix.blocks[i].adjoint();
            }
            return *this;
        }
        BlockDiagonalMatrix& operator*=(const BlockDiagonalMatrixAdjointView<Number>& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] *= rhs._matrix.blocks[i].adjoint();
            }
            return *this;
        }
	
        template<class __matrix__>
        BlockDiagonalMatrix& operator+=(const Eigen::DiagonalWrapper<__matrix__>& rhs) {
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] += rhs.diagonal().segment(this->blocks_begin[i], this->blocks[i].rows()).asDiagonal();
            }
            return *this;
        }
        template<class __matrix__>
        BlockDiagonalMatrix& operator-=(const Eigen::DiagonalWrapper<__matrix__>& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] -= rhs.diagonal().segment(this->blocks_begin[i], this->blocks[i].rows()).asDiagonal();
            }
            return *this;
        }
        template<class __matrix__>
        BlockDiagonalMatrix& operator*=(const Eigen::DiagonalWrapper<__matrix__>& rhs) {
#pragma omp parallel for
            for(int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] *= rhs.diagonal().segment(this->blocks_begin[i], this->blocks[i].rows()).asDiagonal();
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
    BlockDiagonalMatrix<Number> operator+(BlockDiagonalMatrix<Number> lhs, const BlockDiagonalMatrixAdjointView<Number>& rhs) {
        return (lhs += rhs);
    }
    template <class Number>
    BlockDiagonalMatrix<Number> operator-(BlockDiagonalMatrix<Number> lhs, const BlockDiagonalMatrixAdjointView<Number>& rhs) {
        return (lhs -= rhs);
    }
    template <class Number>
    BlockDiagonalMatrix<Number> operator*(BlockDiagonalMatrix<Number> lhs, const BlockDiagonalMatrixAdjointView<Number>& rhs) {
        return (lhs *= rhs);
    }

    template <class Number>
    BlockDiagonalMatrix<Number> operator+(const BlockDiagonalMatrixAdjointView<Number>& lhs, BlockDiagonalMatrix<Number> rhs) {
        return (rhs += lhs);
    }
    template <class Number>
    BlockDiagonalMatrix<Number> operator-(const BlockDiagonalMatrixAdjointView<Number>& lhs, BlockDiagonalMatrix<Number> rhs) {
        return -(rhs -= lhs);
    }
    template <class Number>
    BlockDiagonalMatrix<Number> operator*(const BlockDiagonalMatrixAdjointView<Number>& lhs, BlockDiagonalMatrix<Number> rhs) {
#pragma omp parallel for
        for(int i = 0; i < rhs.blocks.size(); ++i) {
            rhs.blocks[i].applyOnTheLeft(lhs._matrix.blocks[i].adjoint());
        }
        return rhs;
    }

    template <class Number, class __matrix__>
    BlockDiagonalMatrix<Number> operator+(BlockDiagonalMatrix<Number> lhs, const Eigen::DiagonalWrapper<__matrix__>& rhs) {
        return (lhs += rhs);
    }
    template <class Number, class __matrix__>
    BlockDiagonalMatrix<Number> operator-(BlockDiagonalMatrix<Number> lhs, const Eigen::DiagonalWrapper<__matrix__>& rhs) {
        return (lhs -= rhs);
    }
    template <class Number, class __matrix__>
    BlockDiagonalMatrix<Number> operator*(BlockDiagonalMatrix<Number> lhs, const Eigen::DiagonalWrapper<__matrix__>& rhs) {
        return (lhs *= rhs);
    }

    template <class Number, class __matrix__>
    BlockDiagonalMatrix<Number> operator+(const Eigen::DiagonalWrapper<__matrix__>& lhs, BlockDiagonalMatrix<Number> rhs) {
        return (rhs += lhs);
    }
    template <class Number, class __matrix__>
    BlockDiagonalMatrix<Number> operator-(const Eigen::DiagonalWrapper<__matrix__>& lhs, BlockDiagonalMatrix<Number> rhs) {
        return -(rhs -= lhs);
    }
    template <class Number, class __matrix__>
    BlockDiagonalMatrix<Number> operator*(const Eigen::DiagonalWrapper<__matrix__>& lhs, BlockDiagonalMatrix<Number> rhs) {
#pragma omp parallel for
        for(int i = 0; i < rhs.blocks.size(); ++i) {
            rhs.blocks[i].applyOnTheLeft(lhs.diagonal().segment(rhs.blocks_begin[i], rhs.blocks[i].rows()).asDiagonal());
        }
        return rhs;
    }

    template <class EigenMatrixType, class BlockNumber>
    EigenMatrixType operator*(EigenMatrixType basic_matrix, const BlockDiagonalMatrix<BlockNumber>& block_matrix) {
        assert(basic_matrix.cols() == block_matrix.rows());
#pragma omp parallel for
        for(int i = 0; i < block_matrix.blocks.size(); ++i){
            basic_matrix.block(0, block_matrix.blocks_begin[i], 
                basic_matrix.rows(), block_matrix.blocks[i].cols()).applyOnTheRight(block_matrix.blocks[i]);
        }
        return basic_matrix;
    }
    template <class EigenMatrixType, class BlockNumber>
    EigenMatrixType operator*(const BlockDiagonalMatrix<BlockNumber>& block_matrix, EigenMatrixType basic_matrix) {
        assert(block_matrix.cols() == basic_matrix.rows());
#pragma omp parallel for
        for(int i = 0; i < block_matrix.blocks.size(); ++i){
            basic_matrix.block(block_matrix.blocks_begin[i], 0, 
                block_matrix.blocks[i].rows(), basic_matrix.cols()).applyOnTheLeft(block_matrix.blocks[i]);
        }
        return basic_matrix;
    }

    template <class EigenMatrixType, class BlockNumber>
    EigenMatrixType operator*(EigenMatrixType basic_matrix, const BlockDiagonalMatrixAdjointView<BlockNumber>& block_view) {
        assert(basic_matrix.cols() == block_view.rows());
#pragma omp parallel for
        for(int i = 0; i < block_view._matrix.blocks.size(); ++i){
            basic_matrix.block(0, block_view._matrix.blocks_begin[i], 
                basic_matrix.rows(), block_view._matrix.blocks[i].rows()).applyOnTheRight(block_view._matrix.blocks[i].adjoint());
        }
        return basic_matrix;
    }
    template <class EigenMatrixType, class BlockNumber>
    EigenMatrixType operator*(const BlockDiagonalMatrixAdjointView<BlockNumber>& block_view, EigenMatrixType basic_matrix) {
        assert(block_view.cols() == basic_matrix.rows());
#pragma omp parallel for
        for(int i = 0; i < block_view._matrix.blocks.size(); ++i){
            basic_matrix.block(block_view._matrix.blocks_begin[i], 0, 
                block_view._matrix.blocks[i].cols(), basic_matrix.cols()).applyOnTheLeft(block_view._matrix.blocks[i].adjoint());
        }
        return basic_matrix;
    }

}
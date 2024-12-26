#pragma once
#include <Eigen/Dense>
#include <vector>
#include <type_traits>
#include "../UnderlyingFloatingPoint.hpp"

//#define _BLOCKS_USE_OMP

namespace mrock::Utility::Numerics {
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

    template<class _matrix_base>
    struct base_traits {
        using mutable_base = std::remove_cvref_t<_matrix_base>;
        using mutable_vector = std::vector<mutable_base>;
        using as_is_vector = std::conditional_t<std::is_const_v<_matrix_base>, std::add_const_t<mutable_vector>, mutable_vector>;
    };

    template<class source, class target>
    concept not_same_as = !std::is_same_v<target, source>;

    template<class _matrix_base>
    struct BlockDiagonalMatrix {
        using InternalBase = _matrix_base;
        using Scalar = InternalBase::Scalar;
        using ConstructedMatrix = detail::MatrixN<Scalar>;

        typename base_traits<_matrix_base>::as_is_vector blocks;
        std::vector<Eigen::Index> blocks_begin;

        using adjoint_view = decltype(blocks.front().adjoint());
        using eval_result = decltype(blocks.front().eval());

        BlockDiagonalMatrix() = default;
        BlockDiagonalMatrix(const ConstructedMatrix& matrix, const std::vector<HermitianBlock>& block_indizes) 
        {
            blocks.reserve(block_indizes.size());
            blocks_begin.reserve(block_indizes.size());
            for (const auto& block_index : block_indizes)
            {
                blocks.push_back(matrix.block(block_index.position, block_index.position, block_index.size, block_index.size));
                blocks_begin.push_back(block_index.position);
            }
        }
        explicit BlockDiagonalMatrix(const ConstructedMatrix& matrix) 
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
        BlockDiagonalMatrix(typename base_traits<_matrix_base>::mutable_vector&& _blocks, std::vector<Eigen::Index> const& _blocks_begin) 
            : blocks(std::move(_blocks)), blocks_begin(_blocks_begin) {};

        template<not_same_as<_matrix_base> _other_base>
        BlockDiagonalMatrix(BlockDiagonalMatrix<_other_base> const& other)
            : blocks(other.blocks.size()), blocks_begin(other.blocks_begin)
        {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < other.blocks.size(); ++i) {
                this->blocks[i] = other.blocks[i].eval();
            }
        }

        template<not_same_as<_matrix_base> _other_base>
        BlockDiagonalMatrix& operator=(BlockDiagonalMatrix<_other_base> const& other)
        {
            this->blocks.resize(other.blocks.size());
            this->blocks_begin = other.blocks_begin;
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < other.blocks.size(); ++i) {
                this->blocks[i] = other.blocks[i].eval();
            }
            return *this;
        }

        ConstructedMatrix construct_matrix() const 
        {
            ConstructedMatrix ret = ConstructedMatrix::Zero(rows(), cols());
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < blocks.size(); ++i) {
                ret.block(blocks_begin[i], blocks_begin[i], blocks[i].rows(), blocks[i].cols()) = blocks[i];
            }
            return ret;
        }

        BlockDiagonalMatrix<eval_result> eval() const {
            typename base_traits<eval_result>::mutable_vector new_blocks(this->blocks.size());
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < this->blocks.size(); ++i) {
                new_blocks[i] = this->blocks[i].eval();
            }
            return BlockDiagonalMatrix<eval_result>(std::move(new_blocks), blocks_begin);
        }

        BlockDiagonalMatrix<adjoint_view> adjoint() const {
            typename base_traits<adjoint_view>::mutable_vector new_blocks;
            new_blocks.reserve(this->blocks.size());
            for (const auto& block : this->blocks) {
                new_blocks.push_back(block.adjoint());
            }
            return BlockDiagonalMatrix<adjoint_view>(std::move(new_blocks), blocks_begin);
        }
        BlockDiagonalMatrix& adjointInPlace() {
            for (auto& block : this->blocks) {
                block.adjointInPlace();
            }
            return *this;
        }

        inline const InternalBase& block(size_t i) const {
            return blocks[i];
        }
        inline InternalBase& block(size_t i) {
            return blocks[i];
        }

        inline Eigen::Index rows() const {
            Eigen::Index __rows{};
            for (const auto& block : blocks) {
                __rows += block.rows();
            }
            return __rows;
        }
        inline Eigen::Index cols() const {
            Eigen::Index __cols{};
            for (const auto& block : blocks) {
                __cols += block.cols();
            }
            return __cols;
        }

        BlockDiagonalMatrix& applyOnTheLeft(BlockDiagonalMatrix const& other) {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (size_t i = 0U; i < blocks.size(); ++i) {
                this->blocks[i].applyOnTheLeft(other.blocks[i]);
            }
            return *this;
        }

        BlockDiagonalMatrix operator-() const {
            BlockDiagonalMatrix copy(*this);
            for (auto& block : copy.blocks) {
                block *= Scalar{-1.};
            }
            return copy;
        }

        template<class _other_base>
        BlockDiagonalMatrix& operator+=(const BlockDiagonalMatrix<_other_base>& rhs) {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] += rhs.blocks[i];
            }
            return *this;
        }

        template<class _other_base>
        BlockDiagonalMatrix& operator-=(const BlockDiagonalMatrix<_other_base>& rhs) {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] -= rhs.blocks[i];
            }
            return *this;
        }

        template<class _other_base>
        BlockDiagonalMatrix& operator*=(const BlockDiagonalMatrix<_other_base>& rhs) {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] *= rhs.blocks[i];
            }
            return *this;
        }

        template<class __matrix__>
        BlockDiagonalMatrix& operator+=(const Eigen::DiagonalWrapper<__matrix__>& rhs) {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] += rhs.diagonal().segment(this->blocks_begin[i], this->blocks[i].rows()).asDiagonal();
            }
            return *this;
        }

        template<class __matrix__>
        BlockDiagonalMatrix& operator-=(const Eigen::DiagonalWrapper<__matrix__>& rhs) {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] -= rhs.diagonal().segment(this->blocks_begin[i], this->blocks[i].rows()).asDiagonal();
            }
            return *this;
        }

        template<class __matrix__>
        BlockDiagonalMatrix& operator*=(const Eigen::DiagonalWrapper<__matrix__>& rhs) {
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
            for (int i = 0; i < blocks.size(); ++i) {
                this->blocks[i] *= rhs.diagonal().segment(this->blocks_begin[i], this->blocks[i].rows()).asDiagonal();
            }
            return *this;
        }
    };

    namespace detail {
        template<class _base, class _other_base>
        using addition_result = decltype(std::declval<_base>() + std::declval<_other_base>());
        template<class _base, class _other_base>
        using substraction_result = decltype(std::declval<_base>() - std::declval<_other_base>());
        template<class _base, class _other_base>
        using multiplication_result = decltype(std::declval<_base>() * std::declval<_other_base>());
    }

    /*
    *   Basic operations between two BlockDiagonalMatrix objects
    */

    template<class _base, class _other_base>
    BlockDiagonalMatrix<detail::addition_result<_base, _other_base>> operator+(BlockDiagonalMatrix<_base> const& lhs, BlockDiagonalMatrix<_other_base> const& rhs) {
        typename base_traits<detail::addition_result<_base, _other_base>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (size_t i = 0U; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] + rhs.blocks[i]);
        }
        return BlockDiagonalMatrix<detail::addition_result<_base, _other_base>>(std::move(new_blocks), lhs.blocks_begin);
    }

    template<class _base, class _other_base>
    BlockDiagonalMatrix<detail::substraction_result<_base, _other_base>> operator-(BlockDiagonalMatrix<_base> const& lhs, BlockDiagonalMatrix<_other_base> const& rhs) {
        typename base_traits<detail::substraction_result<_base, _other_base>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (size_t i = 0U; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] - rhs.blocks[i]);
        }
        return BlockDiagonalMatrix<detail::substraction_result<_base, _other_base>>(std::move(new_blocks), lhs.blocks_begin);
    }

    template<class _base, class _other_base>
    BlockDiagonalMatrix<detail::multiplication_result<_base, _other_base>> operator*(BlockDiagonalMatrix<_base> const& lhs, BlockDiagonalMatrix<_other_base> const& rhs) {
        typename base_traits<detail::multiplication_result<_base, _other_base>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (size_t i = 0U; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] * rhs.blocks[i]);
        }
        return BlockDiagonalMatrix<detail::multiplication_result<_base, _other_base>>(std::move(new_blocks), lhs.blocks_begin);
    }

    /*
    *   Basic operations between Eigen::Matrix objects and BlockDiagonalMatrix objects
    */

    template <class EigenMatrixType, class _base>
    EigenMatrixType operator*(EigenMatrixType basic_matrix, const BlockDiagonalMatrix<_base>& block_matrix) {
        assert(basic_matrix.cols() == block_matrix.rows());
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < block_matrix.blocks.size(); ++i) {
            basic_matrix.block(0, block_matrix.blocks_begin[i], 
                basic_matrix.rows(), block_matrix.blocks[i].cols()).applyOnTheRight(block_matrix.blocks[i]);
        }
        return basic_matrix;
    }
    template <class EigenMatrixType, class _base>
    EigenMatrixType operator*(const BlockDiagonalMatrix<_base>& block_matrix, EigenMatrixType basic_matrix) {
        assert(block_matrix.cols() == basic_matrix.rows());
#ifdef _BLOCKS_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < block_matrix.blocks.size(); ++i) {
            basic_matrix.block(block_matrix.blocks_begin[i], 0, 
                block_matrix.blocks[i].rows(), basic_matrix.cols()).applyOnTheLeft(block_matrix.blocks[i]);
        }
        return basic_matrix;
    }

    /*
    *   Basic operations between Eigen::DiagonalWrapper and BlockDiagonalMatrix
    */
    namespace detail {
        template<class _matrix>
        using diagonal_segment_type = decltype(std::declval<Eigen::DiagonalWrapper<_matrix>>().diagonal().segment(0, 1).asDiagonal());
    }

    template <class _matrix_base, class __matrix__>
    BlockDiagonalMatrix<detail::addition_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>> operator+(BlockDiagonalMatrix<_matrix_base> const& lhs, Eigen::DiagonalWrapper<__matrix__> const& rhs) {
        typename base_traits<detail::addition_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (int i = 0; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] + rhs.diagonal().segment(lhs.blocks_begin[i], lhs.blocks[i].rows()).asDiagonal());
        }
        return BlockDiagonalMatrix<detail::addition_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>>(std::move(new_blocks), lhs.blocks_begin);
    }

    template <class _matrix_base, class __matrix__>
    BlockDiagonalMatrix<detail::substraction_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>> operator-(BlockDiagonalMatrix<_matrix_base> const& lhs, Eigen::DiagonalWrapper<__matrix__> const& rhs) {
        typename base_traits<detail::substraction_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (int i = 0; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] - rhs.diagonal().segment(lhs.blocks_begin[i], lhs.blocks[i].rows()).asDiagonal());
        }
        return BlockDiagonalMatrix<detail::substraction_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>>(std::move(new_blocks), lhs.blocks_begin);
    }

    template <class _matrix_base, class __matrix__>
    BlockDiagonalMatrix<detail::multiplication_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>> operator*(BlockDiagonalMatrix<_matrix_base> const& lhs, Eigen::DiagonalWrapper<__matrix__> const& rhs) {
        typename base_traits<detail::multiplication_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (int i = 0; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] * rhs.diagonal().segment(lhs.blocks_begin[i], lhs.blocks[i].rows()).asDiagonal());
        }
        return BlockDiagonalMatrix<detail::multiplication_result<_matrix_base, detail::diagonal_segment_type<__matrix__>>>(std::move(new_blocks), lhs.blocks_begin);
    }
}
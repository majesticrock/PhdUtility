#pragma once
#include <Eigen/Dense>
#include <vector>
#include <type_traits>

#include "UnderlyingRealType.hpp"

//#define _BLOCKS_USE_OMP

namespace mrock::iEoM::detail {
    template<class Number> 
    using MatrixN = Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>;
    template<class EigenMatrixType> 
    using RealScalar = UnderlyingRealType_t<typename EigenMatrixType::Scalar>;

    template<class _base, class _other_base>
    using addition_result = decltype(std::declval<_base>() + std::declval<_other_base>());
    template<class _base, class _other_base>
    using substraction_result = decltype(std::declval<_base>() - std::declval<_other_base>());
    template<class _base, class _other_base>
    using multiplication_result = decltype(std::declval<_base>() * std::declval<_other_base>());

    template<class _matrix>
    using diagonal_segment_type = decltype(std::declval<Eigen::DiagonalWrapper<_matrix>>().diagonal().segment(0, 1).asDiagonal());

    /**
     * @brief Describes a Hermitian block within a pivoted matrix.
     *
     * Specifies the starting position and block size for a contiguous
     * Hermitian submatrix.
     */
    struct HermitianBlock {
		Eigen::Index position{};
		Eigen::Index size{};
	};

    /**
     * @brief Stream output helper for HermitianBlock.
     *
     * @param os Output stream.
     * @param block Hermitian block to stringify.
     * @return Reference to the output stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const HermitianBlock& block) {
		os << block.position << "\t" << block.size;
		return os;
	};

    /**
     * @brief Identify contiguous Hermitian blocks in a matrix.
     *
     * @tparam EigenMatrixType Matrix type to inspect.
     * @param matrix Matrix to analyze.
     * @param epsilon Tolerance below which off-diagonal entries are treated as zero.
     * @return Vector of HermitianBlock descriptors for each identified block.
     */
	template<class EigenMatrixType>
	std::vector<HermitianBlock> identify_hermitian_blocks(const EigenMatrixType& matrix, const RealScalar<EigenMatrixType> epsilon = 1e-12) {
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

    /**
     * @brief Traits for block storage type handling.
     *
     * Provides the mutable vector type for a given matrix block type and
     * respects const qualification on the base matrix type.
     */
    template<class _matrix_base>
    struct base_traits {
        using mutable_base = std::remove_cvref_t<_matrix_base>;
        using mutable_vector = std::vector<mutable_base>;
        using as_is_vector = std::conditional_t<std::is_const_v<_matrix_base>, std::add_const_t<mutable_vector>, mutable_vector>;
    };

    /**
     * @brief Concept requiring two types to be different.
     */
    template<class source, class target>
    concept not_same_as = !std::is_same_v<target, source>;

    /**
     * @brief Block diagonal matrix view for a collection of submatrices.
     *
     * This class stores a set of disjoint diagonal blocks and supports
     * block-wise algebraic operations, evaluation, and matrix construction.
     *
     * @tparam _matrix_base Matrix type used for each block.
     */
    template<class _matrix_base>
    struct BlockDiagonalMatrix {
        using InternalBase = _matrix_base;
        using Scalar = InternalBase::Scalar;
        using ConstructedMatrix = MatrixN<Scalar>;

        typename base_traits<_matrix_base>::as_is_vector blocks;
        std::vector<Eigen::Index> blocks_begin;

        using adjoint_view = decltype(blocks.front().adjoint());
        using eval_result = decltype(blocks.front().eval());

        BlockDiagonalMatrix() = default;

        /**
         * @brief Construct from a full matrix and explicit Hermitian block positions.
         *
         * @param matrix Full matrix containing diagonal blocks.
         * @param block_indizes Block descriptors indicating positions and sizes.
         */
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

        /**
         * @brief Construct from a full matrix by automatically identifying Hermitian blocks.
         *
         * @param matrix Full matrix to analyze and partition.
         */
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

        /**
         * @brief Construct directly from a prepared block list and block offsets.
         *
         * @param _blocks Precomputed block matrices.
         * @param _blocks_begin Starting indices for each block.
         */
        BlockDiagonalMatrix(typename base_traits<_matrix_base>::mutable_vector&& _blocks, std::vector<Eigen::Index> const& _blocks_begin) 
            : blocks(std::move(_blocks)), blocks_begin(_blocks_begin) {};

        /**
         * @brief Construct from another BlockDiagonalMatrix with a different block type.
         *
         * @tparam _other_base Block matrix type of the source object.
         * @param other Source block diagonal matrix.
         */
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

        /**
         * @brief Assign from another BlockDiagonalMatrix with a different block type.
         *
         * @tparam _other_base Block matrix type of the source object.
         * @param other Source block diagonal matrix.
         * @return Reference to this object.
         */
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

        /**
         * @brief Reconstruct the full matrix from its diagonal blocks.
         *
         * @return Dense matrix assembled from the contained blocks.
         */
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

        /**
         * @brief Evaluate all block expressions and return a concrete matrix view.
         *
         * @return BlockDiagonalMatrix with evaluated block contents.
         */
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

        /**
         * @brief Return a block diagonal matrix of adjoint blocks.
         *
         * @return BlockDiagonalMatrix containing the adjoint of each block.
         */
        BlockDiagonalMatrix<adjoint_view> adjoint() const {
            typename base_traits<adjoint_view>::mutable_vector new_blocks;
            new_blocks.reserve(this->blocks.size());
            for (const auto& block : this->blocks) {
                new_blocks.push_back(block.adjoint());
            }
            return BlockDiagonalMatrix<adjoint_view>(std::move(new_blocks), blocks_begin);
        }

        /**
         * @brief Replace each block by its adjoint in place.
         *
         * @return Reference to this object.
         */
        BlockDiagonalMatrix& adjointInPlace() {
            for (auto& block : this->blocks) {
                block.adjointInPlace();
            }
            return *this;
        }

        /**
         * @brief Access a block by index.
         *
         * @param i Zero-based block index.
         * @return Const reference to the requested block.
         */
        inline const InternalBase& block(size_t i) const {
            return blocks[i];
        }
        /**
         * @brief Access a block by index.
         *
         * @param i Zero-based block index.
         * @return Reference to the requested block.
         */
        inline InternalBase& block(size_t i) {
            return blocks[i];
        }

        /**
         * @brief Total row dimension of the reconstructed matrix.
         *
         * @return Sum of the row sizes of all blocks.
         */
        inline Eigen::Index rows() const {
            Eigen::Index __rows{};
            for (const auto& block : blocks) {
                __rows += block.rows();
            }
            return __rows;
        }
        /**
         * @brief Total column dimension of the reconstructed matrix.
         *
         * @return Sum of the column sizes of all blocks.
         */
        inline Eigen::Index cols() const {
            Eigen::Index __cols{};
            for (const auto& block : blocks) {
                __cols += block.cols();
            }
            return __cols;
        }

        /**
         * @brief Apply another block diagonal matrix on the left.
         *
         * Each block is transformed by the corresponding block from @p other.
         *
         * @param other Block diagonal matrix to apply.
         * @return Reference to this block diagonal matrix.
         */
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

        /**
         * @brief Add another block diagonal matrix elementwise.
         *
         * @tparam _other_base Block matrix type of the right-hand side.
         * @param rhs Matrix to add.
         * @return Reference to this object.
         */
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

        /**
         * @brief Subtract another block diagonal matrix elementwise.
         *
         * @tparam _other_base Block matrix type of the right-hand side.
         * @param rhs Matrix to subtract.
         * @return Reference to this object.
         */
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

        /**
         * @brief Multiply each corresponding block by another block matrix.
         *
         * @tparam _other_base Block matrix type of the right-hand side.
         * @param rhs Matrix to multiply.
         * @return Reference to this object.
         */
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

        /**
         * @brief Add a diagonal matrix to the block diagonal matrix.
         *
         * Only the diagonal segments overlapping each block are added.
         *
         * @tparam __matrix__ Matrix type wrapped by Eigen::DiagonalWrapper.
         * @param rhs Diagonal matrix wrapper.
         * @return Reference to this object.
         */
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

        /**
         * @brief Subtract a diagonal matrix from the block diagonal matrix.
         *
         * Only the diagonal segments overlapping each block are subtracted.
         *
         * @tparam __matrix__ Matrix type wrapped by Eigen::DiagonalWrapper.
         * @param rhs Diagonal matrix wrapper.
         * @return Reference to this object.
         */
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

        /**
         * @brief Multiply blocks by a diagonal matrix.
         *
         * Only the diagonal segments overlapping each block are used.
         *
         * @tparam __matrix__ Matrix type wrapped by Eigen::DiagonalWrapper.
         * @param rhs Diagonal matrix wrapper.
         * @return Reference to this object.
         */
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

    /*
    *   Basic operations between two BlockDiagonalMatrix objects
    */

    /**
     * @brief Add two block diagonal matrices elementwise.
     *
     * @tparam _base Block type of the left-hand side.
     * @tparam _other_base Block type of the right-hand side.
     * @param lhs Left-hand side block diagonal matrix.
     * @param rhs Right-hand side block diagonal matrix.
     * @return New block diagonal matrix representing the sum.
     */
    template<class _base, class _other_base>
    BlockDiagonalMatrix<addition_result<_base, _other_base>> operator+(BlockDiagonalMatrix<_base> const& lhs, BlockDiagonalMatrix<_other_base> const& rhs) {
        typename base_traits<addition_result<_base, _other_base>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (size_t i = 0U; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] + rhs.blocks[i]);
        }
        return BlockDiagonalMatrix<addition_result<_base, _other_base>>(std::move(new_blocks), lhs.blocks_begin);
    }

    /**
     * @brief Subtract two block diagonal matrices elementwise.
     *
     * @tparam _base Block type of the left-hand side.
     * @tparam _other_base Block type of the right-hand side.
     * @param lhs Left-hand side block diagonal matrix.
     * @param rhs Right-hand side block diagonal matrix.
     * @return New block diagonal matrix representing the difference.
     */
    template<class _base, class _other_base>
    BlockDiagonalMatrix<substraction_result<_base, _other_base>> operator-(BlockDiagonalMatrix<_base> const& lhs, BlockDiagonalMatrix<_other_base> const& rhs) {
        typename base_traits<substraction_result<_base, _other_base>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (size_t i = 0U; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] - rhs.blocks[i]);
        }
        return BlockDiagonalMatrix<substraction_result<_base, _other_base>>(std::move(new_blocks), lhs.blocks_begin);
    }

    /**
     * @brief Multiply two block diagonal matrices elementwise.
     *
     * @tparam _base Block type of the left-hand side.
     * @tparam _other_base Block type of the right-hand side.
     * @param lhs Left-hand side block diagonal matrix.
     * @param rhs Right-hand side block diagonal matrix.
     * @return New block diagonal matrix representing the product.
     */
    template<class _base, class _other_base>
    BlockDiagonalMatrix<multiplication_result<_base, _other_base>> operator*(BlockDiagonalMatrix<_base> const& lhs, BlockDiagonalMatrix<_other_base> const& rhs) {
        typename base_traits<multiplication_result<_base, _other_base>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (size_t i = 0U; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] * rhs.blocks[i]);
        }
        return BlockDiagonalMatrix<multiplication_result<_base, _other_base>>(std::move(new_blocks), lhs.blocks_begin);
    }

    /*
    *   Basic operations between Eigen::Matrix objects and BlockDiagonalMatrix objects
    */

    /**
     * @brief Multiply a dense matrix on the right by a block diagonal matrix.
     *
     * @tparam EigenMatrixType Dense Eigen matrix type.
     * @tparam _base Block type of the block diagonal matrix.
     * @param basic_matrix Dense matrix on the left.
     * @param block_matrix Block diagonal matrix on the right.
     * @return Result of the multiplication.
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
    /**
     * @brief Multiply a block diagonal matrix on the left by a dense matrix.
     *
     * @tparam EigenMatrixType Dense Eigen matrix type.
     * @tparam _base Block type of the block diagonal matrix.
     * @param block_matrix Block diagonal matrix on the left.
     * @param basic_matrix Dense matrix on the right.
     * @return Result of the multiplication.
     */
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
    
    /**
     * @brief Add a diagonal wrapper to a block diagonal matrix.
     *
     * Only diagonal elements that overlap the block structure are added.
     *
     * @tparam _matrix_base Block type of the left-hand side.
     * @tparam __matrix__ Matrix type wrapped by Eigen::DiagonalWrapper.
     * @param lhs Block diagonal matrix.
     * @param rhs Diagonal wrapper to add.
     * @return New block diagonal matrix representing the sum.
     */
    template <class _matrix_base, class __matrix__>
    BlockDiagonalMatrix<addition_result<_matrix_base, diagonal_segment_type<__matrix__>>> operator+(BlockDiagonalMatrix<_matrix_base> const& lhs, Eigen::DiagonalWrapper<__matrix__> const& rhs) {
        typename base_traits<addition_result<_matrix_base, diagonal_segment_type<__matrix__>>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (int i = 0; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] + rhs.diagonal().segment(lhs.blocks_begin[i], lhs.blocks[i].rows()).asDiagonal());
        }
        return BlockDiagonalMatrix<addition_result<_matrix_base, diagonal_segment_type<__matrix__>>>(std::move(new_blocks), lhs.blocks_begin);
    }

    /**
     * @brief Subtract a diagonal wrapper from a block diagonal matrix.
     *
     * Only diagonal elements that overlap the block structure are subtracted.
     *
     * @tparam _matrix_base Block type of the left-hand side.
     * @tparam __matrix__ Matrix type wrapped by Eigen::DiagonalWrapper.
     * @param lhs Block diagonal matrix.
     * @param rhs Diagonal wrapper to subtract.
     * @return New block diagonal matrix representing the difference.
     */
    template <class _matrix_base, class __matrix__>
    BlockDiagonalMatrix<substraction_result<_matrix_base, diagonal_segment_type<__matrix__>>> operator-(BlockDiagonalMatrix<_matrix_base> const& lhs, Eigen::DiagonalWrapper<__matrix__> const& rhs) {
        typename base_traits<substraction_result<_matrix_base, diagonal_segment_type<__matrix__>>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (int i = 0; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] - rhs.diagonal().segment(lhs.blocks_begin[i], lhs.blocks[i].rows()).asDiagonal());
        }
        return BlockDiagonalMatrix<substraction_result<_matrix_base, diagonal_segment_type<__matrix__>>>(std::move(new_blocks), lhs.blocks_begin);
    }

    /**
     * @brief Multiply a block diagonal matrix by a diagonal wrapper.
     *
     * Only diagonal elements that overlap the block structure are used.
     *
     * @tparam _matrix_base Block type of the left-hand side.
     * @tparam __matrix__ Matrix type wrapped by Eigen::DiagonalWrapper.
     * @param lhs Block diagonal matrix.
     * @param rhs Diagonal wrapper to multiply.
     * @return New block diagonal matrix representing the product.
     */
    template <class _matrix_base, class __matrix__>
    BlockDiagonalMatrix<multiplication_result<_matrix_base, diagonal_segment_type<__matrix__>>> operator*(BlockDiagonalMatrix<_matrix_base> const& lhs, Eigen::DiagonalWrapper<__matrix__> const& rhs) {
        typename base_traits<multiplication_result<_matrix_base, diagonal_segment_type<__matrix__>>>::mutable_vector new_blocks;
        new_blocks.reserve(lhs.blocks.size());
        for (int i = 0; i < lhs.blocks.size(); ++i) {
            new_blocks.push_back(lhs.blocks[i] * rhs.diagonal().segment(lhs.blocks_begin[i], lhs.blocks[i].rows()).asDiagonal());
        }
        return BlockDiagonalMatrix<multiplication_result<_matrix_base, diagonal_segment_type<__matrix__>>>(std::move(new_blocks), lhs.blocks_begin);
    }
}
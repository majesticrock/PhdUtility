#ifndef MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_GRAMSCHMIDT_HPP
#define MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_GRAMSCHMIDT_HPP
#include <Eigen/Dense>
#include <limits>
#include <vector>

namespace mrock::iEoM::detail
{
	template <typename DataType>
	class GramSchmidt {
		typedef Eigen::Vector<DataType, Eigen::Dynamic> vector_t;

		/**
		 * @brief Compute the projection of one normalized vector onto another.
		 *
		 * @param from Vector to project.
		 * @param to Unit-length basis vector to project onto.
		 * @return Projection of @p from onto @p to.
		 *
		 * This helper assumes @p to is normalized.
		 */
		inline static vector_t projection(const vector_t& from, const vector_t& to) {
			return from.dot(to) * to;
		};
	public:
		/**
		 * @brief Orthogonalize a vector against an existing orthonormal basis.
		 *
		 * @param vector Vector to be orthogonalized in place.
		 * @param basis Orthonormal basis vectors.
		 * @return Reference to the orthogonalized vector.
		 */
		static vector_t& orthogonalize_single_vector(vector_t& vector, const std::vector<vector_t>& basis) {
			for (const auto& basis_vec : basis)
			{
				vector -= projection(vector, basis_vec);
			}
			return vector;
		};

		/**
		 * @brief Compute an orthonormal basis from the input vectors.
		 *
		 * @param vectors Input vectors to orthogonalize.
		 * @return New vector list containing the orthonormalized vectors.
		 *
		 * The first vector is normalized and subsequent vectors are Gram-Schmidt
		 * orthogonalized against the previous basis vectors.
		 */
		static std::vector<vector_t> compute(const std::vector<vector_t>& vectors) {
			auto buffer = vectors;
			buffer.front().normalize();
			for (std::size_t i = 1U; i < vectors.size(); ++i)
			{
				buffer[i] = vectors[i];
				for (std::size_t j = 0U; j < i; ++j)
				{
					buffer[i] -= projection(buffer[i], buffer[j]);
				}
				buffer[i].normalize();
			}

			return buffer;
		};

		/**
		 * @brief Orthogonalize and normalize the input vectors in place.
		 *
		 * @param vectors Vector list to transform.
		 * @return Reference to the orthogonalized input vector list.
		 */
		static std::vector<vector_t>& compute_and_overwrite(std::vector<vector_t>& vectors) {
			vectors.front().normalize();
			for (std::size_t i = 1U; i < vectors.size(); ++i)
			{
				for (std::size_t j = 0U; j < i; ++j)
				{
					vectors[i] -= projection(vectors[i], vectors[j]);
				}
				vectors[i].normalize();
			}
			return vectors;
		};
	};
}
#endif  // MROCK_IEOM_INCLUDE_MROCK_IEOM_DETAIL_GRAMSCHMIDT_HPP

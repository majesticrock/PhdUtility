/**
 * @file MomentumList.hpp
 * @brief Defines the MomentumList class for handling a list of Momentum objects.
 */

#pragma once
#include "Momentum.hpp"
#include <mrock/utility/VectorWrapper.hpp>
#include <algorithm>

namespace mrock::symbolic_operators {

	/**
	 * @class MomentumList
	 * @brief A wrapper class for a vector of Momentum objects with additional functionalities.
	 */
	struct MomentumList : public mrock::utility::VectorWrapper<Momentum>
	{
	private:
		using _parent = mrock::utility::VectorWrapper<Momentum>;

	public:
		/**
		 * @brief Serializes the MomentumList object, required for boost support.
		 * @tparam Archive The type of the archive.
		 * @param ar The archive object.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->_vector;
		}

		/**
		 * @brief Default constructor.
		 */
		MomentumList();

		/**
		 * @brief Constructs a MomentumList with a single Momentum object.
		 * @param momentum The Momentum object to initialize the list with.
		 */
		explicit MomentumList(const Momentum& momentum);

		/**
		 * @brief Constructs a MomentumList with two Momentum objects.
		 * @param first The first Momentum object.
		 * @param second The second Momentum object.
		 */
		MomentumList(const Momentum& first, const Momentum& second);

		/**
		 * @brief Constructs a MomentumList with an initializer list of Momentum objects.
		 * @param init The initializer list of Momentum objects.
		 */
		MomentumList(std::initializer_list<Momentum> init);

		/**
		 * @brief Constructs a MomentumList with an initializer list of characters.
		 * @param init The initializer list of characters.
		 */
		MomentumList(std::initializer_list<char> init);

		/**
		 * @brief Multiplies each Momentum object in the list by a given factor.
		 * @param rhs The factor to multiply by.
		 * @return A reference to the modified MomentumList.
		 */
		inline MomentumList& operator*=(const int rhs);

		/**
		 * @brief Multiplies each Momentum object in the list by a given factor.
		 * @param factor The factor to multiply by.
		 */
		inline void multiply_by(int factor);

		/**
		 * @brief Flips the momentum of each Momentum object in the list.
		 */
		inline void flip_momentum();

		/**
		 * @brief Sorts the Momentum objects in the list.
		 */
		inline void sort();

		/**
		 * @brief Replaces occurrences of a specific MomentumSymbol name with a given Momentum object.
		 * @param replaceWhat The MomentumSymbol name to replace.
		 * @param replaceWith The Momentum object to replace with.
		 */
		void replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith);

		/**
		 * @brief Removes Momentum objects with zero value from the list.
		 */
		void remove_zeros();

		/**
		 * @brief Flips the momentum of a single Momentum object identified by its MomentumSymbol name.
		 * @param momentum The MomentumSymbol name of the Momentum object to flip.
		 */
		void flip_single(const MomentumSymbol::name_type momentum);
	};

	/**
	 * @brief Outputs the MomentumList to an output stream.
	 * @param os The output stream.
	 * @param momenta The MomentumList to output.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const MomentumList& momenta);

	// Inline definitions
	MomentumList& MomentumList::operator*=(const int rhs) {
		for (auto& mom : _vector) {
			mom *= rhs;
		}
		return *this;
	}
	void MomentumList::multiply_by(int factor) {
		(*this) *= factor;
	}
	void MomentumList::flip_momentum() {
		(*this) *= -1;
	}
	void MomentumList::sort() {
		std::sort(this->begin(), this->end());
	}
} // namespace mrock::symbolic_operators
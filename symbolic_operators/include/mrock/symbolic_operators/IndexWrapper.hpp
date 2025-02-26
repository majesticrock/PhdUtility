/**
 * @file IndexWrapper.hpp
 * @brief Defines the Index enum and the IndexWrapper class for handling indizes.
 */

#pragma once
#include <iostream>
#include <mrock/utility/VectorWrapper.hpp>
#include <string>
#include <map>

namespace mrock::symbolic_operators {

	/**
	 * @typedef index_base
	 * @brief Defines the base type for the Index enum as unsigned char.
	 */
	typedef unsigned char index_base;

	/**
	 * @enum Index
	 * @brief Enumeration representing various symbolic indices.
	 * 
	 * The indices include:
	 * - SpinUp: Represents spin up (0).
	 * - SpinDown: Represents spin down (1).
	 * - Sigma: Represents sigma (2).
	 * - SigmaPrime: Represents sigma prime (3).
	 * - GeneralSpin_S: Represents general spin S (4).
	 * - GeneralSpin_SPrime: Represents general spin S prime (5).
	 * - TypeA: Represents type A (6).
	 * - TypeB: Represents type B (7).
	 * - TypeC: Represents type C (8).
	 * - char_a: Represents the lowercase ASCII character 'a' (97).
	 * - UndefinedIndex: Represents an undefined index (254).
	 * - NoIndex: Represents no index (255).
	 * 
	 * - Not explicitly represented, but defined implementation-wise are any lower-case ASCII characters.char_to_index.
	 *	 These can be obtained via static_cast<Index>(char) or using the char_to_index() function
	 */
	enum class Index : index_base { 
		SpinUp = 0, 
		SpinDown, 
		Sigma, 
		SigmaPrime, 
		GeneralSpin_S, 
		GeneralSpin_SPrime, 
		TypeA, 
		TypeB, 
		TypeC,
		char_a = 97, 
		UndefinedIndex = 254, 
		NoIndex = 255 
	};

	/**
	 * @brief Converts a character to an Index.
	 * 
	 * This function is designed for lower-case letters but works for any ASCII representable character.
	 * 
	 * @param c The character to convert.
	 * @return The corresponding Index value.
	 */
	constexpr Index char_to_index(unsigned char c) {
		return static_cast<Index>(c);
	}

	/**
	 * @var string_to_index
	 * @brief A map that associates string representations with their corresponding Index values.
	 */
	inline const std::map<std::string, Index> string_to_index = {
		{"up", Index::SpinUp}, 
		{"down", Index::SpinDown}, 
		{"sigma", Index::Sigma}, 
		{"sigma'", Index::SigmaPrime}, 
		{"S", Index::GeneralSpin_S}, 
		{"S'", Index::GeneralSpin_SPrime}, 
		{"A", Index::TypeA}, 
		{"B", Index::TypeB}, 
		{"C", Index::TypeC},
		{"a", Index::char_a},
		{"b", char_to_index('b')},
		{"c", char_to_index('c')},
		{"d", char_to_index('d')},
		{"e", char_to_index('e')},
		{"f", char_to_index('f')},
		{"g", char_to_index('g')},
		{"h", char_to_index('h')},
		{"i", char_to_index('i')},
		{"j", char_to_index('j')},
		{"k", char_to_index('k')},
		{"l", char_to_index('l')},
		{"m", char_to_index('m')},
		{"n", char_to_index('n')},
		{"o", char_to_index('o')},
		{"p", char_to_index('p')},
		{"q", char_to_index('q')},
		{"r", char_to_index('r')},
		{"s", char_to_index('s')},
		{"t", char_to_index('t')},
		{"u", char_to_index('u')},
		{"v", char_to_index('v')},
		{"w", char_to_index('w')},
		{"x", char_to_index('x')},
		{"y", char_to_index('y')},
		{"z", char_to_index('z')}
	};

	/**
	 * @brief Checks if the given index represents a variable (mutable).
	 * 'Mutable' means that it is associated with a sum or similar.
	 * An example is sigma; it is commonly summed over as a representation of spins.
	 * Then expressions like delta_{sigma,up} can be evaluated to be one if sigma=up.
	 * An Index like SpinUp is set to be non-mutable. This allows us to evaluate delta_{up,down}=0.
	 * 
	 * An index is considered mutable if it is between 2 and 5 (inclusive).
	 * 
	 * @param idx The index to check.
	 * @return True if the index is mutable, false otherwise.
	 */
	constexpr bool is_mutable(const Index idx) {
		return (static_cast<index_base>(idx) > 1 && static_cast<index_base>(idx) < 6);
	}

	/**
	 * @brief Overloads the stream insertion operator for the Index enum.
	 * 
	 * @param os The output stream.
	 * @param index The Index value to insert into the stream.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const Index index);

	/**
	 * @struct IndexWrapper
	 * @brief A wrapper for a vector of Index values.
	 * 
	 * This struct provides serialization support and comparison operators.
	 */
	struct IndexWrapper {
		std::vector<Index> indizes; ///< The vector of Index values.

		/**
		 * @brief Serializes the IndexWrapper, required for boost support.
		 * @tparam Archive The type of the archive.
		 * @param ar The archive to serialize to.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->indizes;
		}

		/**
		 * @brief Default constructor.
		 */
		IndexWrapper() = default;

		/**
		 * @brief Constructs an IndexWrapper with a single Index value.
		 * 
		 * @param _spin The Index value to initialize with.
		 */
		IndexWrapper(Index _spin) : indizes(1U, _spin) {};

		/**
		 * @brief Constructs an IndexWrapper with a vector of Index values.
		 * 
		 * @param _indizes The vector of Index values to initialize with.
		 */
		IndexWrapper(const std::vector<Index>& _indizes)
			: indizes(_indizes) {};

		/**
		 * @brief Constructs an IndexWrapper with a vector of Index values (move semantics).
		 * 
		 * @param _indizes The vector of Index values to initialize with.
		 */
		IndexWrapper(std::vector<Index>&& _indizes)
			: indizes(std::move(_indizes)) {};

		MROCK_VECTOR_WRAPPER_FILL_MEMBERS(Index, indizes);

		/**
		 * @brief Compares two IndexWrapper objects.
		 * 
		 * @param rhs The other IndexWrapper to compare with.
		 * @return The result of the comparison.
		 */
		inline auto operator<=>(const IndexWrapper& rhs) const = default;
	};

	/**
	 * @brief Overloads the stream insertion operator for the IndexWrapper struct.
	 * 
	 * @param os The output stream.
	 * @param indizes The IndexWrapper to insert into the stream.
	 * @return The output stream.
	 */
	std::ostream& operator<<(std::ostream& os, const IndexWrapper& indizes);
} // namespace mrock::symbolic_operators
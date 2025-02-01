#pragma once
#include <iostream>
#include <mrock/utility/VectorWrapper.hpp>
#include <string>
#include <map>

namespace mrock::symbolic_operators {
	typedef unsigned char index_base;
	enum class Index : index_base { SpinUp = 0, SpinDown, Sigma, SigmaPrime, GeneralSpin_S, GeneralSpin_SPrime, 
	TypeA, TypeB, TypeC,
	char_a = 97, // lowercase ASCII characters are represented by their byte value, but we list only the first one explicitly
	UndefinedIndex = 254, NoIndex = 255 };
	inline const std::map<std::string, Index> string_to_index = {
		{"up", Index::SpinUp}, {"down", Index::SpinDown}, {"sigma", Index::Sigma}, {"sigma'", Index::SigmaPrime}, 
		{"S", Index::GeneralSpin_S}, {"S'", Index::GeneralSpin_SPrime}, 
		{"A", Index::TypeA}, {"B", Index::TypeB}, {"C", Index::TypeC},
		{"a", Index::char_a}
	};
	// Turns a character into an Index - is designed for lower-case letters, 
	// but works in principle for any ASCII representable character.
	constexpr Index char_to_index(unsigned char c) {
		return static_cast<Index>(c);
	}
	// Returns true if the index represents a variable and false otherwise
	// Example: If the index is SpinUp it is fixed, i.e., non-mutable
	constexpr bool is_mutable(const Index idx) {
		return (static_cast<index_base>(idx) > 1 && static_cast<index_base>(idx) < 6);
	}

	std::ostream& operator<<(std::ostream& os, const Index index);

	struct IndexWrapper {
		std::vector<Index> indizes;

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& this->indizes;
		}

		IndexWrapper() = default;
		IndexWrapper(Index _spin) : indizes(1U, _spin) {};
		IndexWrapper(const std::vector<Index>& _indizes)
			: indizes(_indizes) {};
		IndexWrapper(std::vector<Index>&& _indizes)
			: indizes(std::move(_indizes)) {};

		VECTOR_WRAPPER_FILL_MEMBERS(Index, indizes);

		inline auto operator<=>(const IndexWrapper& rhs) const = default;
	};

	std::ostream& operator<<(std::ostream& os, const IndexWrapper& indizes);
}
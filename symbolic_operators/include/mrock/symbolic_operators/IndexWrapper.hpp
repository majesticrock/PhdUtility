#pragma once
#include <iostream>
#include <mrock/utility/VectorWrapper.hpp>
#include <string>
#include <map>

namespace mrock::symbolic_operators {
	enum class Index { SpinUp = 0, SpinDown, Sigma, SigmaPrime, GeneralSpin_S, GeneralSpin_SPrime, BosonA, BosonB, UndefinedIndex, NoIndex = 256 };
	inline const std::map<std::string, Index> string_to_index = {
		{"up", Index::SpinUp}, {"down", Index::SpinDown}, {"sigma", Index::Sigma}, {"sigma'", Index::SigmaPrime}, 
		{"S", Index::GeneralSpin_S}, {"S'", Index::GeneralSpin_SPrime}, {"A", Index::BosonA}, {"B", Index::BosonB}
	};
	// Returns true if the index represents a variable and false otherwise
	// Example: If the index is SpinUp it is fixed, i.e., non-mutable
	constexpr bool is_mutable(const Index idx) {
		return (static_cast<int>(idx) > 1 && static_cast<int>(idx) < 6);
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
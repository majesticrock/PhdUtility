#pragma once
#include <iostream>
#include <string>

namespace mrock::symbolic_operators {
    struct MomentumSymbol {
		// Represents essentially nothing but a char, but does not define arithmetic operations
		struct name_type { 
			char _n{};

			template<class Archive>
			void serialize(Archive& ar, const unsigned int version) {
				ar& _n;
			}

			constexpr name_type() = default;
			constexpr name_type(char n) noexcept : _n{n} {};
			constexpr auto operator<=>(const name_type&) const = default;
			constexpr auto operator<=>(const char other) const { return _n <=> other; }

			explicit constexpr operator char() const noexcept { return _n; }
		};

		int factor{};
		name_type name{};

		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& factor;
			ar& name;
		}
		
		constexpr MomentumSymbol() = default;
		constexpr MomentumSymbol(int _factor, char _name) : factor{_factor}, name(_name) {};
        constexpr MomentumSymbol(int _factor, name_type _name) : factor{_factor}, name{_name} {};
		constexpr auto operator<=>(const MomentumSymbol&) const = default;
	};

	inline std::ostream& operator<<(std::ostream& os, const MomentumSymbol::name_type name) { return (os << name._n); }
	inline std::istream& operator>>(std::istream& is, MomentumSymbol::name_type& name) { return (is >> name._n); }

	inline std::string operator+(const std::string& str, const MomentumSymbol::name_type sym) {
		return (str + static_cast<char>(sym));
	}
	inline std::string operator+(const MomentumSymbol::name_type sym, const std::string& str) {
		return (static_cast<char>(sym) + str);
	}
}
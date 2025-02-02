/**
 * @file MomentumSymbol.hpp
 * @brief Defines the MomentumSymbol structure and related operators for symbolic operations.
 */

#pragma once
#include <iostream>
#include <string>

namespace mrock::symbolic_operators {

	/**
	 * @struct MomentumSymbol
	 * @brief Represents a symbolic momentum with a factor and a name.
	 */
	struct MomentumSymbol {
		
		/**
		 * @struct name_type
		 * @brief Represents a name as a single character with comparison and serialization capabilities, 
		 * 		but without arithmetic operations (it does not make sense to add or multiply names).
		 */
		struct name_type { 
			char _n{}; ///< The character representing the name.

			/**
			 * @brief Serializes the name_type object, required for boost support.
			 * @tparam Archive The type of the archive.
			 * @param ar The archive to serialize to.
			 * @param version The version of the serialization.
			 */
			template<class Archive>
			void serialize(Archive& ar, const unsigned int version) {
				ar& _n;
			}

			constexpr name_type() = default;
			/**
			 * @brief Constructs a name_type with a given character.
			 * @param n The character to initialize the name with.
			 */
			constexpr name_type(char n) noexcept : _n{n} {};

			/**
			 * @brief Compares two name_type objects.
			 * @param other The other name_type object to compare with.
			 * @return The result of the comparison.
			 */
			constexpr auto operator<=>(const name_type&) const = default;

			/**
			 * @brief Compares the name_type object with a character.
			 * @param other The character to compare with.
			 * @return The result of the comparison.
			 */
			constexpr auto operator<=>(const char other) const { return _n <=> other; }

			/**
			 * @brief Converts the name_type object to a character.
			 * @return The character representation of the name_type.
			 */
			explicit constexpr operator char() const noexcept { return _n; }
		};

		int factor{}; ///< The factor associated with the momentum.
		name_type name{}; ///< The name associated with the momentum.

		/**
		 * @brief Serializes the MomentumSymbol object, required for boost support.
		 * @tparam Archive The type of the archive.
		 * @param ar The archive to serialize to.
		 * @param version The version of the serialization.
		 */
		template<class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar& factor;
			ar& name;
		}
		
		constexpr MomentumSymbol() = default;

		/**
		 * @brief Constructs a MomentumSymbol with a given factor and name.
		 * @param _factor The factor to initialize the momentum with.
		 * @param _name The name to initialize the momentum with.
		 */
		constexpr MomentumSymbol(int _factor, char _name) : factor{_factor}, name(_name) {};

		/**
		 * @brief Constructs a MomentumSymbol with a given factor and name_type.
		 * @param _factor The factor to initialize the momentum with.
		 * @param _name The name_type to initialize the momentum with.
		 */
		constexpr MomentumSymbol(int _factor, name_type _name) : factor{_factor}, name{_name} {};

		/**
		 * @brief Compares two MomentumSymbol objects.
		 * @param other The other MomentumSymbol object to compare with.
		 * @return The result of the comparison.
		 */
		constexpr auto operator<=>(const MomentumSymbol&) const = default;
	};

	/**
	 * @brief Outputs the name_type to an output stream.
	 * @param os The output stream.
	 * @param name The name_type to output.
	 * @return The output stream.
	 */
	inline std::ostream& operator<<(std::ostream& os, const MomentumSymbol::name_type name) { return (os << name._n); }

	/**
	 * @brief Inputs a name_type from an input stream.
	 * @param is The input stream.
	 * @param name The name_type to input.
	 * @return The input stream.
	 */
	inline std::istream& operator>>(std::istream& is, MomentumSymbol::name_type& name) { return (is >> name._n); }

	/**
	 * @brief Concatenates a string and a name_type.
	 * @param str The string to concatenate.
	 * @param sym The name_type to concatenate.
	 * @return The concatenated string.
	 */
	inline std::string operator+(const std::string& str, const MomentumSymbol::name_type sym) {
		return (str + static_cast<char>(sym));
	}

	/**
	 * @brief Concatenates a name_type and a string.
	 * @param sym The name_type to concatenate.
	 * @param str The string to concatenate.
	 * @return The concatenated string.
	 */
	inline std::string operator+(const MomentumSymbol::name_type sym, const std::string& str) {
		return (static_cast<char>(sym) + str);
	}
} // namespace mrock::symbolic_operators
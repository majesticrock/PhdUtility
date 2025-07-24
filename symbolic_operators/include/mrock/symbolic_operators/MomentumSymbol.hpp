/**
 * @file MomentumSymbol.hpp
 * @brief Defines the MomentumSymbol structure and related operators for symbolic operations.
 */

#ifndef _MROCK_SYM_OP_MOMENTUM_SYMBOL_HPP_
#define _MROCK_SYM_OP_MOMENTUM_SYMBOL_HPP_

#include <iostream>
#include <string>
#include <concepts>

#ifndef MROCK_TEX_VECTOR
/**
 * @def MROCK_TEX_VECTOR
 * @brief A macro that allows users to define their desired vector style during compilation.
 *
 * This macro defines the prefix used for vector notation in string formatting, 
 * allowing customization of how vectors are represented. By default, it is set 
 * to "\\mathbf", but users can override this definition during compilation to 
 * use a different style (e.g., "\\vec", "\\boldsymbol", etc.).
 *
 * Example usage:
 * @code
 * #define MROCK_TEX_VECTOR "\\vec"
 * #include "MomentumSymbol.hpp"
 * @endcode
 */
#define MROCK_TEX_VECTOR "\\mathbf"
#endif

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
 	 * @concept MrockStringAddable
 	 * @brief Specifies that a type can be added to a std::string and result in a std::string.
 	 * @tparam T The type being checked.
 	 */
	template<typename T>
	concept MrockStringAddable = requires(std::string s, T x) {
	    { s + x } -> std::convertible_to<std::string>;
	};

	/**
 	 * @concept MrockValidForVectorWrap
 	 * @brief Specifies that a type can be wrapped with vector notation, excluding MomentumSymbol::name_type.
 	 * @tparam T The type being checked.
 	 */
	template<typename T>
	concept MrockValidForVectorWrap = MrockStringAddable<T> && (!std::same_as<MomentumSymbol::name_type, T>);

	/**
 	 * @concept MrockHasToString
 	 * @brief Specifies that a type has a valid std::to_string implementation.
 	 * @tparam T The type being checked.
 	 */
	template<typename T>
	concept MrockHasToString = requires (T t) {
	    { std::to_string(t) } -> std::convertible_to<std::string>;
	};

	/**
	 * @brief Wraps an object in vector notation using MROCK_TEX_VECTOR prefix.
	 *        This overload is used for types that support concatenation with std::string directly.
	 * @tparam T A type satisfying MrockValidForVectorWrap concept.
	 * @param x The object to wrap in vector notation.
	 * @return A string representing the wrapped object in vector notation.
	 */
	template <MrockValidForVectorWrap T>
	inline std::string _vector_wrap(const T& x)
	{
	    return std::string(MROCK_TEX_VECTOR "{") + x + std::string("}");
	}

	/**
     * @brief Wraps an object in vector notation using MROCK_TEX_VECTOR prefix,
     *        leveraging std::to_string for conversion if supported by the type.
     *        This overload excludes types satisfying MrockValidForVectorWrap concept but supports MrockHasToString concept.
     * @tparam T A type satisfying MrockHasToString, but not MrockValidForVectorWrap.
     * @param x The object to wrap in vector notation.
     * @return A string representing the wrapped object in vector notation.
     */
	template <typename T>
	inline std::string _vector_wrap(const T& x) requires MrockHasToString<T> && (!MrockValidForVectorWrap<T>)
	{
		return std::string(MROCK_TEX_VECTOR "{") + std::to_string(x) + std::string("}");
	}

	/**
     * @brief Wraps a MomentumSymbol::name_type object in vector notation using MROCK_TEX_VECTOR prefix.
     * @param x The name_type object to wrap in vector notation.
     * @return A string representing the wrapped name_type in vector notation.
     */
	inline std::string _vector_wrap(const MomentumSymbol::name_type& x)
	{
    	return std::string(MROCK_TEX_VECTOR "{") + static_cast<char>(x) + "}";
	}

	/**
     * @brief Outputs a MomentumSymbol::name_type object to an output stream, formatted with vector notation.
     * @param os The output stream to write to.
     * @param name The MomentumSymbol::name_type object to output.
     * @return The modified output stream after writing the formatted name_type object.
     */
	inline std::ostream& operator<<(std::ostream& os, const MomentumSymbol::name_type name) { return (os << _vector_wrap(name)); }

	/**
     * @brief Reads a character into a MomentumSymbol::name_type from an input stream.
     *        Updates the internal character of the name_type based on input data.
     * @param is The input stream to read from.
     * @param name The MomentumSymbol::name_type object being updated with input data.
     * @return The modified input stream after reading into the name_type object.
     */
	inline std::istream& operator>>(std::istream& is, MomentumSymbol::name_type& name) { return (is >> name._n); }

	/**
     * @brief Concatenates a std::string and a MomentumSymbol::name_type, formatting the latter with vector notation before concatenation.
     * @param str The std::string to concatenate with the symbol's formatted representation.
     * @param sym The MomentumSymbol::name_type whose formatted representation is concatenated with `str`.
     * @return A new std::string resulting from concatenating `str` and `sym`.
     */
	inline std::string operator+(const std::string& str, const MomentumSymbol::name_type sym) {
		return (str + _vector_wrap(sym));
	}

	/**
     * @brief Concatenates a MomentumSymbol::name_type and a std::string, formatting the former with vector notation before concatenation.
     * @param sym The MomentumSymbol::name_type whose formatted representation is concatenated with `str`.
     * @param str The std::string to concatenate after the symbol's formatted representation.
     * @return A new std::string resulting from concatenating `sym` and `str`.
     */
	inline std::string operator+(const MomentumSymbol::name_type sym, const std::string& str) {
		return (_vector_wrap(sym) + str);
	}
} // namespace mrock::symbolic_operators

#endif // _MROCK_SYM_OP_MOMENTUM_SYMBOL_HPP_
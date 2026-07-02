/**
 * @file string_helper.hpp
 * @brief Small string utility helpers used by the symbolic operators library.
 *
 * This header provides generic utilities for splitting strings, removing
 * escape characters, finding unescaped symbols, and extracting comma-separated
 * elements enclosed by delimiters.
 */
#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <type_traits>
#include <algorithm>

namespace mrock::symbolic_operators {
/**
 * @brief Splits a string into tokens using a single-character delimiter.
 *
 * @tparam StringType A string-like type (e.g. `std::string`).
 * @tparam InputStringStreamType Stream type used for tokenization. Defaults
 *        to a basic_istringstream matching `StringType`.
 * @param str The input string to split.
 * @param delimiter The character used to split the string.
 * @return A vector of tokens of type `StringType`.
 */
	template<class StringType, class InputStringStreamType = std::basic_istringstream<typename StringType::value_type, typename StringType::traits_type, typename StringType::allocator_type>>
	std::vector<StringType> split(const StringType &str, std::add_const_t<typename StringType::value_type> delimiter) {
		std::vector<StringType> tokens;
		StringType token;
		InputStringStreamType tokenStream(str);
		while (std::getline(tokenStream, token, delimiter)) {
			tokens.push_back(token);
		}
		return tokens;
	}
	
/**
 * @brief Removes escape characters from a string in-place.
 *
 * Sequences of double-escape (e.g. `\\\\`) are reduced to a single
 * escape character, while an escape followed by any other character removes
 * the escape and keeps the character. The function returns the number of
 * characters removed.
 *
 * @tparam StringType A string-like type.
 * @tparam ForwardIt Iterator type for the string (defaults to `StringType::iterator`).
 * @param input The string to modify in-place.
 * @param escape The escape character to remove (defaults to '\\').
 * @return The number of characters removed from `input`.
 */
	template<class StringType, class ForwardIt = typename StringType::iterator>
	size_t remove_escape_characters(StringType& input, std::add_const_t<typename StringType::value_type> escape = '\\') {
		assert(!(input.size() == 1U && input.back() == escape) && "The last character of the string must not be an escape character!");
		assert(!(input.size() > 1U && input.back() == escape && *(input.end() - 2) != escape) && "The last character of the string must not be an escape character!");
		
		ForwardIt first = std::find(input.begin(), input.end(), escape);
	    if (first == input.end()) return 0U;

	    for (ForwardIt i = first; i != input.end(); ++i) {
	        if (*i != escape) {
			*first++ = std::move(*i);
		}
		else if ((i + 1) != input.end() && *(i + 1) == escape) {
			// Situation [esc][esc] should become [esc]
			*first++ = std::move(*i++);
		}
	}
	    size_t n_removed = input.end() - first;
		input.erase(first, input.end());
		return n_removed;
	}

/**
 * @brief Finds the first occurrence of a symbol that is not escaped.
 *
 * This behaves like `StringType::find(symbol, start)` but skips any match
 * where the symbol is immediately preceded by the escape character.
 *
 * @tparam StringType A string-like type.
 * @param input The input string to search.
 * @param symbol The character to find.
 * @param start The starting index to search from.
 * @param escape The escape character (defaults to '\\').
 * @return The position of the unescaped symbol, or `StringType::npos` if none found.
 */
	template<class StringType>
	size_t find_skip_escaped(const StringType& input, std::add_const_t<typename StringType::value_type> symbol, size_t start = 0U, 
		std::add_const_t<typename StringType::value_type> escape = '\\') 
	{
		while (start < input.size()) {
			size_t pos = input.find(symbol, start);
			if(pos == StringType::npos) return StringType::npos;
			
			if (pos > 0 && input[pos - 1] == escape) {
				start = pos + 1;
			}
			else {
				return pos;
			}
		}

		return StringType::npos;
	}

/**
 * @brief Extracts comma-separated elements enclosed by delimiters from a string.
 *
 * For example, given `"...{a,b,c}..."` this returns `["a","b","c"]`.
 * The function respects escape characters for the delimiters.
 *
 * @tparam StringType A string-like type.
 * @param input The input string containing the delimited list.
 * @param left_delimiter The opening delimiter (default '{').
 * @param right_delimiter The closing delimiter (default '}').
 * @param escape The escape character used to skip delimiters (default '\\').
 * @return A vector of extracted elements as `StringType` values.
 */
	template<class StringType>
	std::vector<StringType> extract_elements(const StringType& input, std::add_const_t<typename StringType::value_type> left_delimiter = '{', 
		std::add_const_t<typename StringType::value_type> right_delimiter = '}', std::add_const_t<typename StringType::value_type> escape = '\\') 
	{
		std::vector<StringType> elements;
	
		// Find the positions of the braces
		size_t startPos = find_skip_escaped(input, left_delimiter , 0U, escape);
		size_t endPos   = find_skip_escaped(input, right_delimiter, 0U, escape);
	
		// Check if both braces are found
		if (startPos != StringType::npos && endPos != StringType::npos && startPos < endPos) {
			// Extract the substring within the braces
			StringType substring = input.substr(startPos + 1, endPos - startPos - 1);
			elements = split(substring, typename StringType::value_type(','));
		}
	
		return elements;
	}
}
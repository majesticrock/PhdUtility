#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <type_traits>
#include <algorithm>

namespace mrock::utility {
    // Splits the string at delimiter
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

	// Works just as std::string::find but skips the found character, if it is preceeded by an escape character (e.g. \ )
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

	// Extracts elements from a list {x1,x2,x3...} within some string
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
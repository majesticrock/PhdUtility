#pragma once

#include <string>
#include <vector>
#include <sstream>

namespace mrock::utility {
    // Splits the string at delimiter
	inline std::vector<std::string> split(const std::string &str, char delimiter) {
		std::vector<std::string> tokens;
		std::string token;
		std::istringstream tokenStream(str);
		while (std::getline(tokenStream, token, delimiter)) {
			tokens.push_back(token);
		}
		return tokens;
	}
	
	// Extracts elements from a list {x1,x2,x3...} within some string
	inline std::vector<std::string> extract_elements(const std::string &input, char left_delimiter = '{', char right_delimiter = '}') {
		std::vector<std::string> elements;
	
		// Find the positions of the braces
		size_t startPos = input.find(left_delimiter);
		size_t endPos = input.find(right_delimiter);
	
		// Check if both braces are found
		if (startPos != std::string::npos && endPos != std::string::npos && startPos < endPos) {
			// Extract the substring within the braces
			std::string substring = input.substr(startPos + 1, endPos - startPos - 1);
			elements = split(substring, ',');
		}
	
		return elements;
	}
}
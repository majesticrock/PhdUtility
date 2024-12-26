#ifndef _InputFileReader_hpp_
#define _InputFileReader_hpp_

#include <fstream>
#include <vector>
#include <string>

/**
* InputFileReader
* This is a general class to handle input files.
* Each line contains a parameter name and the content of the parameter seperated by a space.
* The content can be a bool, int, double or string. Or a list of them seperated by spaces.
* The order of the parameters is irrelevant
*
* Example file:
*
* method Classic
* M 1024
* ak_dist Box
* N_Box 100
* J_dist Const
* mode Mode1
* t_p 0 0.1 0.4
*
* Example code:
* InputFileReader input("input.config");
* string method = getString("method");
* cout << method << endl // Output: "Classic"
* int M = getInt("M");
* cout << M << endl // Output: 1024
* etc.
**/
namespace mrock::Utility {
	class InputFileReader {
	private:
		// a list of the parameter names
		std::vector<std::string> names;
		// a list of the parameter content
		// this strings have to be decrypted when a getX() method is called
		std::vector<std::string> contents;
		// tracking which parameters are used
		std::vector<bool> used;

		// get the index of the parameter name
		// return -1 if not found
		inline int find(std::string s)
		{
			int pos = -1;
			bool done = false;
			for (unsigned int i = 0U; i < names.size() && !done; ++i) {
				if (names[i] == s) {
					pos = i;
					done = true;
				}
			}
			return pos;
		};

	public:
		// Read an input file with the name "fileName"
		InputFileReader(std::string fileName)
		{
			std::ifstream f;

			f.open(fileName, std::ios::in);

			std::string s;
			while (std::getline(f, s)) {
				if (s[0] == '#') continue;
				size_t pos = s.find(" ");
				std::string name = s.substr(0, pos);
				std::string content = s.substr(pos + 1);

				if (content.back() == 13) {
					content.erase(prev(content.end()));
				}

				names.push_back(name);
				contents.push_back(content);
				used.push_back(false);
			}
		};

		// checks wether the parameter name is specified
		inline bool is(std::string name)
		{
			int i = find(name);
			if (i != -1) {
				used[i] = true;
				return true;
			}
			else {
				return false;
			}
		};
		// get the content as boolean
		inline bool getBool(std::string name)
		{
			int i = find(name);
			if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
			used[i] = true;
			if (contents[i] == "true") {
				return true;
			}
			else if (contents[i] == "false") {
				return false;
			}
			else {
				throw std::invalid_argument("The Parameter " + name + " is not bool; " + contents[i]);
			}
		};
		// get the content as integer
		inline int getInt(std::string name)
		{
			int i = find(name);
			if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
			used[i] = true;
			return std::stoi(contents[i]);
		};
		// get the content as list of integers
		inline std::vector<int> getIntList(std::string name)
		{
			int i = find(name);
			if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
			used[i] = true;

			std::string s = contents[i];
			std::vector<int> result;
			size_t pos = 0;
			bool done = false;
			while (!done) {
				result.push_back(std::stoi(s, &pos));
				if (s.size() > pos) {
					s = s.substr(pos + 1);
				}
				else {
					done = true;
				}
			}
			return result;
		};
		// get the content as double
		inline double getDouble(std::string name)
		{
			int i = find(name);
			if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
			used[i] = true;
			return std::stod(contents[i]);
		};
		// get the content as list of doubles
		inline std::vector<double> getDoubleList(std::string name)
		{
			int i = find(name);
			if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
			used[i] = true;

			std::string s = contents[i];
			std::vector<double> result;
			size_t pos = 0;
			bool done = false;
			while (!done) {
				result.push_back(std::stod(s, &pos));
				if (s.size() > pos) {
					s = s.substr(pos + 1);
				}
				else {
					done = true;
				}
			}
			return result;
		};
		// get the content as string
		inline std::string getString(std::string name)
		{
			int i = find(name);
			if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
			used[i] = true;
			return contents[i];
		};
		// get the content as list of strings
		inline std::vector<std::string> getStringList(std::string name)
		{
			int i = find(name);
			if (i == -1) throw std::invalid_argument("Could not find parameter " + name);
			used[i] = true;
			std::string s = contents[i];
			std::vector<std::string> result;
			bool done = false;
			while (!done) {
				size_t pos2 = s.find(" ");
				if (pos2 == std::string::npos) {
					std::string content = s;
					result.push_back(content);
					done = true;
				}
				else {
					std::string content = s.substr(0, pos2);
					s = s.substr(pos2 + 1);
					result.push_back(content);
				}
			}
			return result;
		};

		// check wether all parameters are used
		// a parameter is used when a is() method is called
		inline bool allUsed()
		{
			bool result = true;
			for (unsigned int i = 0; i < used.size(); i++) {
				if (used[i] == false) result = false;
			}
			return result;
		};
		// list the names of all unused parameters
		inline std::vector<std::string> listNotUsed()
		{
			std::vector<std::string> result;
			for (unsigned int i = 0; i < used.size(); i++) {
				if (!used[i]) {
					result.push_back(names[i]);
				}
			}
			return result;
		};
	};
}
#endif
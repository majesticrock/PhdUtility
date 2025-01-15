#include <string>
#include <iostream>
#include "../include/mrock/utility/StringUtility.hpp"

using namespace mrock::utility;

int main() {
    std::string test = "\\{\\\\a\\\\\\bc\\{{def,ghi,xyz\\}}";

    std::vector<std::string> elements = extract_elements(test);
    std::vector<std::string> expected_elements = {"def", "ghi", "xyz\\}"};
    std::cout << "Found " << elements.size() << " elements. Expected 3." << std::endl;

    for (size_t i = 0U; i < elements.size(); ++i) {
        std::cout << i << ":\t" << elements[i] << "\t" << expected_elements[i] << std::endl;
        if(elements[i] != expected_elements[i]) return 1;
    }
    if (elements.size() != 3U) return 1;

    size_t removed = remove_escape_characters(test);
    std::string expected_string = "{\\a\\bc{{def,ghi,xyz}}";
    std::cout << "Removed " << removed << " characters. Expected 6." << std::endl;
    std::cout << "Resultung string is:     " << test << std::endl;
    std::cout << "Expected:                " << expected_string << std::endl;
    if (removed != 6U) return 1;
    if (test != expected_string) return 1;

    test = "V_\\\\mathrm\\{CUT\\}";
    removed = remove_escape_characters(test);
    expected_string = "V_\\mathrm{CUT}";
    std::cout << "Removed " << removed << " characters. Expected 3." << std::endl;
    std::cout << "Resultung string is:     " << test << std::endl;
    std::cout << "Expected:                " << expected_string << std::endl;
    if (removed != 3U) return 1;
    if (test != expected_string) return 1;

    return 0;
}
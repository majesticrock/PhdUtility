#include <mrock/utility/better_to_string.hpp>

#include <iostream>
#include <string>
#include <utility>

#define FAIL_IF_NOT(cond)                                                      \
    if (!(cond)) {                                                             \
        std::cerr << "FAILED: " #cond << " at line " << __LINE__ << std::endl; \
        return 1;                                                              \
    }

int main() {
    using mrock::utility::better_to_string;
    using mrock::utility::better_to_string_fixed;

    {
        std::cout << "Test 1: integer conversion\n";

        std::string result = better_to_string(42);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 42\n\n";

        FAIL_IF_NOT(result == "42");
    }

    {
        std::cout << "Test 2: negative integer conversion\n";

        std::string result = better_to_string(-12345);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: -12345\n\n";

        FAIL_IF_NOT(result == "-12345");
    }

    {
        std::cout << "Test 3: unsigned integer conversion\n";

        std::string result = better_to_string(123456789u);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 123456789\n\n";

        FAIL_IF_NOT(result == "123456789");
    }

    {
        std::cout << "Test 4: integer conversion with base 16\n";

        std::string result = better_to_string(255, 16);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: ff\n\n";

        FAIL_IF_NOT(result == "ff");
    }

    {
        std::cout << "Test 5: integer conversion with base 2\n";

        std::string result = better_to_string(10, 2);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 1010\n\n";

        FAIL_IF_NOT(result == "1010");
    }

    {
        std::cout << "Test 6: fixed floating-point conversion\n";

        std::string result = better_to_string_fixed(3.141592653589793, 3);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 3.142\n\n";

        FAIL_IF_NOT(result == "3.142");
    }

    {
        std::cout << "Test 7: fixed floating-point conversion with zero digits\n";

        std::string result = better_to_string_fixed(3.75, 0);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: 4\n\n";

        FAIL_IF_NOT(result == "4");
    }

    {
        std::cout << "Test 8: negative fixed floating-point conversion\n";

        std::string result = better_to_string_fixed(-12.3456, 2);

        std::cout << "Result:   " << result << '\n';
        std::cout << "Expected: -12.35\n\n";

        FAIL_IF_NOT(result == "-12.35");
    }

    {
        std::cout << "Test 9: buffer-too-small error handling\n";

        // The internal buffer has only 32 chars.
        // Requesting 100 fixed digits should fail and return an error message.
        std::string result = better_to_string_fixed(1.2345, 100);

        std::cout << "Resulting error message: " << result << "\n\n";

        FAIL_IF_NOT(!result.empty());
        FAIL_IF_NOT(result != "1.2345");
    }

    std::cout << "All better_to_string tests passed.\n";

    return 0;
}
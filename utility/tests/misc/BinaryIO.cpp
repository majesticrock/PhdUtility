#include <mrock/utility/BinaryIO.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#define FAIL_IF_NOT(cond)                                                      \
    if (!(cond)) {                                                             \
        std::cerr << "FAILED: " #cond << " at line " << __LINE__ << std::endl; \
        return 1;                                                              \
    }

int main() {
    using namespace mrock::utility::BinaryIO;

    const std::string variable_file = "test_binaryio_variable.bin";
    const std::string vector_file = "test_binaryio_vector.bin";
    const std::string serialized_vector_file = "test_binaryio_serialized_vector.bin";
    const std::string empty_vector_file = "test_binaryio_empty_vector.bin";

    {
        std::cout << "Test 1: writeVariable / read_to_variable\n";

        const int original = 123456;
        int loaded = 0;

        {
            auto writer = create_writer(variable_file);
            FAIL_IF_NOT(writer.is_open());

            writeVariable(original, writer);
            writer.close();
        }

        {
            auto reader = create_reader(variable_file);
            FAIL_IF_NOT(reader.is_open());

            read_to_variable(loaded, reader);
            reader.close();
        }

        std::cout << "Original value: " << original << '\n';
        std::cout << "Loaded value:   " << loaded << "\n\n";

        FAIL_IF_NOT(loaded == original);
    }

    {
        std::cout << "Test 2: writeVector / readToVector\n";

        const std::vector<double> original{1.25, -2.5, 3.75, 42.0};
        std::vector<double> loaded(original.size());

        {
            auto writer = create_writer(vector_file);
            FAIL_IF_NOT(writer.is_open());

            writeVector(original, writer);
            writer.close();
        }

        {
            auto reader = create_reader(vector_file);
            FAIL_IF_NOT(reader.is_open());

            readToVector(loaded, reader);
            reader.close();
        }

        std::cout << "Original vector size: " << original.size() << '\n';
        std::cout << "Loaded vector size:   " << loaded.size() << "\n\n";

        FAIL_IF_NOT(loaded.size() == original.size());

        for (std::size_t i = 0; i < original.size(); ++i) {
            FAIL_IF_NOT(loaded[i] == original[i]);
        }
    }

    {
        std::cout << "Test 3: serializeVector / readSerializedVector\n";

        const std::vector<int> original{1, 2, 3, 4, 5, -10, 12345};
        std::vector<int> loaded;

        {
            auto writer = serializeVector(original, serialized_vector_file);
            FAIL_IF_NOT(writer.is_open());
            writer.close();
        }

        {
            auto reader = readSerializedVector(loaded, serialized_vector_file);
            FAIL_IF_NOT(reader.is_open());
            reader.close();
        }

        std::cout << "Original vector size: " << original.size() << '\n';
        std::cout << "Loaded vector size:   " << loaded.size() << "\n\n";

        FAIL_IF_NOT(loaded.size() == original.size());

        for (std::size_t i = 0; i < original.size(); ++i) {
            FAIL_IF_NOT(loaded[i] == original[i]);
        }
    }

    {
        std::cout << "Test 4: serializeVector / readSerializedVector with empty vector\n";

        const std::vector<float> original{};
        std::vector<float> loaded{1.0f, 2.0f, 3.0f};

        {
            auto writer = serializeVector(original, empty_vector_file);
            FAIL_IF_NOT(writer.is_open());
            writer.close();
        }

        {
            auto reader = readSerializedVector(loaded, empty_vector_file);
            FAIL_IF_NOT(reader.is_open());
            reader.close();
        }

        std::cout << "Original vector size: " << original.size() << '\n';
        std::cout << "Loaded vector size:   " << loaded.size() << "\n\n";

        FAIL_IF_NOT(loaded.empty());
    }

    {
        std::cout << "Test 5: opening non-existing file for reading should fail\n";

        std::vector<int> loaded;
        auto reader = readSerializedVector(loaded, "this_file_should_not_exist.bin");

        FAIL_IF_NOT(!reader.is_open() || !reader.good());

        std::cout << "Non-existing file test passed.\n\n";
    }

    std::remove(variable_file.c_str());
    std::remove(vector_file.c_str());
    std::remove(serialized_vector_file.c_str());
    std::remove(empty_vector_file.c_str());

    std::cout << "All BinaryIO tests passed.\n";

    return 0;
}
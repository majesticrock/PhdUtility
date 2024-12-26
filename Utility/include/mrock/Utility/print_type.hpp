#pragma once
#include <iostream>

namespace mrock::Utility::Debug {
    template<class T>
    void print_type() {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }
}
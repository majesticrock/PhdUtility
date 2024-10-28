#pragma once
#include <iostream>

namespace Utility::Debug {
    template<class T>
    void print_type() {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }
}
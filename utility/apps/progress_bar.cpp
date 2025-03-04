/**
 * @file progress_bar.cpp
 * @brief Example for the usage of the progress bar utility
 * 
 * Compile with g++ -fopenmp progress_bar.cpp -o progress_bar.o
 * Run with ./progress_bar.o
 */

#include "../include/mrock/utility/progress_bar.hpp"
#include <chrono>
#include <thread>
#include <omp.h>
#include <cstdlib>

// from https://en.cppreference.com/w/cpp/numeric/random/rand
unsigned bounded_rand(unsigned range)
{
    for (unsigned x, r;;)
        if (x = rand(), r = x % range, x - r <= -range)
            return r;
}


int main() {
    for (int i = 0; i < 100; ++i) {
        mrock::utility::progress_bar(i / 100.f);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    std::cout << std::endl;

    std::vector<float> progresses(16);
#pragma omp parallel for
    for (int n = 0; n < 16; ++n) {
        for (int i = 0; i < 100; ++i) {
            progresses[n] = i / 100.f;
            if (omp_get_thread_num() == 0) {
                mrock::utility::multi_progress_bar(16, progresses);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(60U + bounded_rand(120U)));
        }
    }
    std::cout << std::endl;
}
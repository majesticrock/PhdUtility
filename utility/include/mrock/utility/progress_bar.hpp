#pragma once

#include <iostream>
#include <vector>

namespace mrock::utility {
    inline void multi_progress_bar(int n, std::vector<float> const& progresses) {
        const int bar_length = 128 / n;
        std::cout << "[";

        float total_progress{};
        for (int i = 0; i < n; ++i) {
            total_progress += progresses[i] / n;
            const int position = bar_length * progresses[i];
            for (int j = 0; j < bar_length; ++j) {
                if (j < position) {
                    std::cout << "|";
                }
                else if (j == position) {
                    std::cout << ">";
                }
                else {
                    std::cout << " ";
                }
            }
        }
        std::cout << "]" << int(100.f * total_progress) << "%\r";
        std::cout.flush();
    }

    inline void progress_bar(float progress) {
        constexpr int bar_length = 50;
        const int position = bar_length * progress;
        for (int j = 0; j < bar_length; ++j) {
            if (j < position) {
                std::cout << "|";
            }
            else if (j == position) {
                std::cout << ">";
            }
            else {
                std::cout << " ";
            }
        }
        std::cout << "]" << int(100.f * progress) << "%\r";
        std::cout.flush();
    }
}
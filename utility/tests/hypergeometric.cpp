#include "../include/mrock/utility/Numerics/hypergeometric_2F1.hpp"

#include <iostream>
#include <iomanip>
#include <complex>

template<class T>
bool check(T numeric, T compare) {
    std::cout << "We computed: " << numeric << "   WolframAlpha computed: " << compare;
    std::cout << "    || Error: " << std::abs(compare - numeric) << std::endl;

    return std::abs(compare - numeric) < 1e-7;
}

int main() {
    using namespace mrock::utility::Numerics;

    std::cout << std::scientific << std::setprecision(12);

    { // Real arguments
        std::complex<double> result = hypergeometric_2F1(0.2, 0.8, 0.3, 0.4);
        std::complex<double> wolfram_result = 1.3245222276888917546110595458;
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5, 0.5, 1.0, 0.4);
        wolfram_result = 1.131603977657727948718251611721;
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(2.0, 1.5, 1.0, 0.9);
        wolfram_result = 458.5302607244150031398395639427;
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5, 1.0, 1.0, 2.4);
        wolfram_result = {0.0, -0.84515425472851657750961832736594 };
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5L, 0.5L, 1.0L, 2.4L, 1e-10L);
        wolfram_result = {0.7352710849571768692031125,  -0.7939111082485524424133790314399};
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5L, 0.5L, 1.0L, 0.9999L, 1e-12L);
        wolfram_result = 3.814364242073590871011732087;
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5L, 0.5L, 1.0L, 1.0001L, 1e-12L);
        wolfram_result = {3.814205358821811209128288881588, -0.99997500140615235122620107017173};
        if (!check(result, wolfram_result)) return 1;

        // This is too close to the branch cut |z| and the routine fails.
        //result = hypergeometric_2F1(0.5L, 0.5L, 1.0L, 1.00000001L, 1e-12L);
        //wolfram_result = {6.7460271763725098367934, -0.99999999750000001406};
        //if (!check(result, wolfram_result)) return 1;
    }

    { // Real arguments
        std::complex<double> result = hypergeometric_2F1(0.2, 0.8, 0.3, std::complex<double>{0.4, 0.4});
        std::complex<double> wolfram_result = {1.11909047589423097919477333, 0.379370335626998712098792957913};
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5, 0.5, 1.0, std::complex<double>{-0.4, -0.4});
        wolfram_result = {0.90702832764897288538482118, -0.0661847419908417524511};
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(2.0, 1.5, 1.0, std::complex<double>{0., 0.9});
        wolfram_result = {-0.33009753622503674186024186, 0.40480048598360060470562156699066};
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5, 1.0, 1.0, std::complex<double>{2.4, 2.4});
        wolfram_result = {0.2987981970454975064162992, 0.5202187163268035594653567986};
        if (!check(result, wolfram_result)) return 1;

        result = hypergeometric_2F1(0.5L, 0.5L, 1.0L, std::complex<long double>{-2.4L, 2.4L}, 1e-10L);
        wolfram_result = {0.67116253785715921597323, 0.1224626024182782632274};
        if (!check(result, wolfram_result)) return 1;
    }

    return 0;
}
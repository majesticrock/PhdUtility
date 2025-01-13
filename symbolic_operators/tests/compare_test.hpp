#pragma once
#include "../include/mrock/symbolic_operators/Term.hpp"
#include "../include/mrock/symbolic_operators/Wick.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>
#include <filesystem>
#include <string>

namespace sym_op_test {
    using namespace mrock::symbolic_operators;

    struct SymOpTest {
        const std::string COMPARE_DIR;

        SymOpTest(const std::string _compare_dir) : COMPARE_DIR(_compare_dir) {};

        template <class TestClass>
        bool perform_comparison(const TestClass& A, const TestClass& B) 
        {
            if (A.size() != B.size()) {
                return false;
            }
            for (size_t i = 0U; i < A.size(); ++i) {
                if (A[i] != B[i]) {
                    std::cerr << "No match for: \nA: " << A[i] << "\nB: " << B[i] << std::endl;
                    return false;
                }
            }
            return true;
        }

        // Returns true if the test was successful
        template <class TestClass>
        bool load_and_test(const std::string& name, const TestClass& computed) 
        {
            std::ifstream ifs(COMPARE_DIR + name, std::ios::binary);
    	    if(ifs.good()) {
                TestClass compare;
    	    	boost::archive::binary_iarchive ia(ifs);
    	    	ia >> compare;
    	    	ifs.close();
                if(!perform_comparison(computed, compare)) {
                    std::cerr << "Failed comparison at " << name << " !" << std::endl;
                    return false;
                }
    	    }
    	    else {
    	    	throw std::runtime_error("Inputstream for " + COMPARE_DIR + name + " is bad!");
    	    }
            return true;
        }

        template <class TestClass>
        void save_as_comparison(const std::string& name, const TestClass& correct) 
        {
            std::filesystem::create_directories(COMPARE_DIR);
            std::ofstream ofs(COMPARE_DIR + name, std::ios::binary);
            boost::archive::binary_oarchive oa(ofs);
            oa << correct;
            ofs.close();
        }

        inline int perform_test(const std::vector<Term>& H, const std::vector<Term>& base_term, 
            const std::vector<WickOperatorTemplate>& templates, const std::vector<std::unique_ptr<WickSymmetry>>& symmetries,
            const bool is_baseline)
        {
            const std::string file_names[3] = { "first_commutation.bin", "second_commutation.bin", "wicks.bin" };
            /*
            *
            *   Computing and testing of [H, B]
            * 
            */
            std::vector<Term> first_commutation;
            commutator(first_commutation, H, base_term);
            clean_up(first_commutation);
            // This means the test has been passed. We generate a new comparison file
            if (is_baseline) save_as_comparison(file_names[0], first_commutation);
            if(!load_and_test(file_names[0], first_commutation)) return 1;

            /*
            *
            *   Computing and testing of [B^+, [H, B]]
            * 
            */
            std::vector<Term> second_commutation;
            std::vector<Term> base_dagger = base_term;
            rename_momenta(base_dagger, 'k', 'l');
            hermitian_conjugate(base_dagger);
            commutator(second_commutation, base_dagger, first_commutation);
            clean_up(second_commutation);
            // This means the test has been passed. We generate a new comparison file
            if (is_baseline) save_as_comparison(file_names[1], second_commutation);
            if(!load_and_test(file_names[1], second_commutation)) return 1;

            /*
            *
            *   Computing and testing of Wick's theorem on <[B^+, [H, B]]>
            * 
            */
            WickTermCollector wicks;
            wicks_theorem(second_commutation, templates, wicks);
            clear_etas(wicks);
            clean_wicks(wicks, symmetries);
            // This means the test has been passed. We generate a new comparison file
            if (is_baseline) save_as_comparison(file_names[2], wicks);
            if(!load_and_test(file_names[2], wicks)) return 1;

            return 0;
        }
    };
}
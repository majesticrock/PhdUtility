#include <mrock/symbolic_operators/ExpectationValues>

#include <vector>

using namespace mrock::symbolic_operators;

int main() {
    const std::vector<WickOperatorTemplate> number_template({
         /* Template for number operators */
         WickOperatorTemplate{
             { Num_Comparison }, Momentum(), OperatorType::Number
         }
    });
    const std::vector<WickOperatorTemplate> number_and_sc_template({
         /* Template for number operators */
         WickOperatorTemplate{
             { Num_Comparison }, Momentum(), OperatorType::Number
         },
         /* Template for pair creation/annihilation operators */
         WickOperatorTemplate{
             { SC_Comparison }, Momentum(), OperatorType::SC
         }
    });

    const Momentum base_k = Momentum('k');

    // Setup a few operators to be used later
    const Operator c_k = Operator{base_k, Index::SpinUp, false};
    const Operator c_minus_k = Operator{-base_k, Index::SpinDown, false};
    const Operator c_k_dagger = Operator{base_k, Index::SpinUp, true};
    const Operator c_minus_k_dagger = Operator{-base_k, Index::SpinDown, true};


    const TermCollector pair_creation({ Term(
        1, std::vector<Operator>({
            c_k_dagger, c_minus_k_dagger.with_momentum('p')
        })
    )});

    WickTermCollector only_num_wick = wicks_theorem(pair_creation, number_template);
    only_num_wick.clean_up();
    WickTermCollector num_and_sc_wick = wicks_theorem(pair_creation, number_and_sc_template);
    num_and_sc_wick.clean_up();
    
    std::cout << "If only number expectation values are finite\n\\begin{align*}\n\t"
        "\\langle" << pair_creation << "\\rangle = " << only_num_wick << "\\end{align*}"
        << std::endl;

    std::cout << "If BCS-like expectation values are also finite\n\\begin{align*}\n\t"
        "\\langle" << pair_creation << "\\rangle = " << num_and_sc_wick << "\\end{align*}"
        << std::endl;

    return 0;
}
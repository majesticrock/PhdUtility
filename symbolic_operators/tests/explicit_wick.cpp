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

    const WickTermCollector expected_only_num_wick;
    WickTerm expected_sc_term;
    expected_sc_term.multiplicity = 1;
    expected_sc_term.delta_momenta.push_back({ Momentum('p'), Momentum('k') });
    expected_sc_term.operators.emplace_back(OperatorType::SC, true, Momentum('k'));
    const WickTermCollector expected_num_and_sc_wick({ expected_sc_term });

    WickTerm doubled_sc_term = expected_sc_term;
    doubled_sc_term.multiplicity = 2;
    WickTermCollector duplicate_wicks({ expected_sc_term, doubled_sc_term });
    duplicate_wicks.combine_duplicates();

    TermCollector duplicate_terms({ Term(1, std::vector<Operator>{ c_k }),
                                    Term(2, std::vector<Operator>{ c_k }) });
    duplicate_terms.combine_duplicates();
    
    std::cout << "If only number expectation values are finite\n\\begin{align*}\n\t"
        "\\langle" << pair_creation << "\\rangle = " << only_num_wick << "\\end{align*}"
        << std::endl;

    std::cout << "If BCS-like expectation values are also finite\n\\begin{align*}\n\t"
        "\\langle" << pair_creation << "\\rangle = " << num_and_sc_wick << "\\end{align*}"
        << std::endl;

    if (only_num_wick.terms != expected_only_num_wick.terms ||
        num_and_sc_wick.terms != expected_num_and_sc_wick.terms || duplicate_wicks.size() != 1U ||
        duplicate_wicks.front().multiplicity != 3 || duplicate_terms.size() != 1U ||
        duplicate_terms.front().multiplicity != 3) {
        return 1;
    }
    return 0;
}
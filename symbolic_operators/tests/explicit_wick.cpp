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
    const Momentum transfer_q('q');
    const std::vector<WickOperatorTemplate> zero_and_finite_momentum_templates({
        WickOperatorTemplate{{ Num_Comparison }, Momentum(), OperatorType::Number},
        WickOperatorTemplate{{ SC_Comparison }, Momentum(), OperatorType::SC},
        WickOperatorTemplate{{ Num_Comparison }, -transfer_q, OperatorType::CDW},
        WickOperatorTemplate{{ SC_Comparison }, transfer_q, OperatorType::Eta}
    });

    const Momentum base_k = Momentum('k');
    const Momentum base_p = Momentum('p');

    // Setup a few operators to be used later
    const Operator c_k = Operator{base_k, Index::SpinUp, false};
    const Operator c_minus_k = Operator{-base_k, Index::SpinDown, false};
    const Operator c_k_dagger = Operator{base_k, Index::SpinUp, true};
    const Operator c_minus_k_dagger = Operator{-base_k, Index::SpinDown, true};
    const Operator c_p = Operator{base_p, Index::SpinUp, false};
    const Operator c_minus_p = Operator{-base_p, Index::SpinDown, false};


    const TermCollector pair_creation({ Term(
        1, std::vector<Operator>({
            c_k_dagger, c_minus_k_dagger.with_momentum('p')
        })
    )});

    const TermCollector pair_annihilation({ Term(1, std::vector<Operator>{ c_minus_k, c_k }) });
    const TermCollector number_term({ Term(1, std::vector<Operator>{ c_k_dagger, c_p }) });
    const TermCollector reversed_number_term({ Term(1, std::vector<Operator>{ c_p, c_k_dagger }) });
    const TermCollector cdw_term({ Term(1, std::vector<Operator>{ c_k_dagger, c_p }) });
    const TermCollector eta_term({ Term(1, std::vector<Operator>{ c_minus_p, c_k }) });

    WickTermCollector only_num_wick = wicks_theorem(pair_creation, number_template);
    only_num_wick.clean_up();
    WickTermCollector num_and_sc_wick = wicks_theorem(pair_creation, number_and_sc_template);
    num_and_sc_wick.clean_up();

    WickTermCollector cc_wick = wicks_theorem(pair_annihilation, number_and_sc_template);
    cc_wick.clean_up();
    WickTermCollector c_dagger_c_wick = wicks_theorem(number_term, number_and_sc_template);
    c_dagger_c_wick.clean_up();
    WickTermCollector c_c_dagger_wick = wicks_theorem(reversed_number_term, zero_and_finite_momentum_templates);
    c_c_dagger_wick.clean_up();
    WickTermCollector cdw_wick = wicks_theorem(cdw_term, zero_and_finite_momentum_templates);
    cdw_wick.clean_up();
    WickTermCollector eta_wick = wicks_theorem(eta_term, zero_and_finite_momentum_templates);
    eta_wick.clean_up();

    const WickTermCollector expected_only_num_wick;
    WickTerm expected_sc_term;
    expected_sc_term.multiplicity = 1;
    expected_sc_term.delta_momenta.push_back({ Momentum('p'), Momentum('k') });
    expected_sc_term.operators.emplace_back(OperatorType::SC, true, Momentum('k'));
    const WickTermCollector expected_num_and_sc_wick({ expected_sc_term });

    WickTerm expected_cc_term;
    expected_cc_term.operators.emplace_back(OperatorType::SC, false, Momentum('k'));
    const WickTermCollector expected_cc_wick({ expected_cc_term });
    WickTerm expected_c_dagger_c_term;
    expected_c_dagger_c_term.operators.emplace_back(OperatorType::Number, false, Momentum('k'), Index::SpinUp);
    expected_c_dagger_c_term.delta_momenta.push_back({ base_p, base_k });
    const WickTermCollector expected_c_dagger_c_wick({ expected_c_dagger_c_term });
    WickTerm expected_c_c_dagger_term;
    expected_c_c_dagger_term.multiplicity = 1;
    expected_c_c_dagger_term.operators.emplace_back(OperatorType::HoleNumber, false, Momentum('k'), Index::SpinUp);
    expected_c_c_dagger_term.delta_momenta.push_back({ base_p, base_k });
    WickTerm expected_c_c_dagger_cdw_term;
    expected_c_c_dagger_cdw_term.multiplicity = -1;
    expected_c_c_dagger_cdw_term.operators.emplace_back(OperatorType::CDW, false, base_p, Index::SpinUp);
    expected_c_c_dagger_cdw_term.delta_momenta.push_back({ base_p, base_k + transfer_q });
    const WickTermCollector expected_c_c_dagger_wick({ expected_c_c_dagger_cdw_term, expected_c_c_dagger_term });
    WickTerm expected_cdw_term;
    expected_cdw_term.operators.emplace_back(OperatorType::CDW, false, base_p, Index::SpinUp);
    expected_cdw_term.delta_momenta.push_back({ base_p, base_k + transfer_q });
    WickTermCollector expected_cdw_wick({ expected_cdw_term, expected_c_dagger_c_term });
    WickTerm expected_zero_momentum_eta_term;
    expected_zero_momentum_eta_term.operators.emplace_back(OperatorType::SC, false, base_k);
    expected_zero_momentum_eta_term.delta_momenta.push_back({ base_p, base_k });
    WickTerm expected_eta_term;
    expected_eta_term.operators.emplace_back(OperatorType::Eta, false, base_p - transfer_q);
    expected_eta_term.delta_momenta.push_back({ base_p, base_k + transfer_q });
    WickTermCollector expected_eta_wick({ expected_eta_term, expected_zero_momentum_eta_term });

    WickTerm doubled_sc_term = expected_sc_term;
    doubled_sc_term.multiplicity = 2;
    WickTermCollector duplicate_wicks({ expected_sc_term, doubled_sc_term });
    duplicate_wicks.combine_duplicates();

    TermCollector duplicate_terms({ Term(1, std::vector<Operator>{ c_k }),
                                    Term(2, std::vector<Operator>{ c_k }) });
    duplicate_terms.combine_duplicates();
    
    std::cout << "% If only number expectation values are finite\n\\begin{align*}\n\t"
        "\\langle" << pair_creation << "\\rangle = " << only_num_wick << "\\end{align*}"
        << std::endl;

    std::cout << "% If BCS-like expectation values are also finite\n\\begin{align*}\n\t"
        "\\langle" << pair_creation << "\\rangle = " << num_and_sc_wick << "\\end{align*}"
        << std::endl;

    std::cout << "% Wick's theorem for $cc$\n\\begin{align*}\n\t"
        "\\langle" << pair_annihilation << "\\rangle = " << cc_wick << "\\end{align*}"
        << std::endl;
    std::cout << "% Wick's theorem for $c^{\\dagger} c$\n\\begin{align*}\n\t"
        "\\langle" << number_term << "\\rangle = " << c_dagger_c_wick << "\\end{align*}"
        << std::endl;
    std::cout << "% Wick's theorem for $c c^{\\dagger}$\n\\begin{align*}\n\t"
        "\\langle" << reversed_number_term << "\\rangle = " << c_c_dagger_wick << "\\end{align*}"
        << std::endl;
    std::cout << "% Wick's theorem with zero- and finite-momentum number-like expectations\n\\begin{align*}\n\t"
        "\\langle" << cdw_term << "\\rangle = " << cdw_wick << "\\end{align*}"
        << std::endl;
    std::cout << "% Wick's theorem with zero- and finite-momentum pair expectations\n\\begin{align*}\n\t"
        "\\langle" << eta_term << "\\rangle = " << eta_wick << "\\end{align*}"
        << std::endl;

    if (only_num_wick.terms != expected_only_num_wick.terms ||
        num_and_sc_wick.terms != expected_num_and_sc_wick.terms || duplicate_wicks.size() != 1U ||
        cc_wick.terms != expected_cc_wick.terms ||
        c_dagger_c_wick.terms != expected_c_dagger_c_wick.terms ||
        c_c_dagger_wick.size() != 2U || c_c_dagger_wick.terms != expected_c_c_dagger_wick.terms ||
        cdw_wick.size() != 2U || eta_wick.size() != 2U ||
        cdw_wick.terms != expected_cdw_wick.terms || eta_wick.terms != expected_eta_wick.terms ||
        duplicate_wicks.front().multiplicity != 3 || duplicate_terms.size() != 1U ||
        duplicate_terms.front().multiplicity != 3) {
        return 1;
    }
    return 0;
}
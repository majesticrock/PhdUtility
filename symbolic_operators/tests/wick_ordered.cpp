#include <mrock/symbolic_operators/experimental/WickOrderedCollector.hpp>

#include <mrock/symbolic_operators/Term.hpp>

using namespace mrock::symbolic_operators;

int main(int /*argc*/, char** /*argv*/) {
    const Operator c_q_up     = Operator(Momentum('q'), Index::SpinUp,   false);
    const Operator c_q_down   = Operator(Momentum('q'), Index::SpinDown, false);
    const Operator c_q_sigma  = Operator(Momentum('q'), Index::Sigma,    false);

    const Term H_kin(1, 
        Coefficient::RealInversionSymmetric("\\varepsilon", MomentumList(Momentum('q'))),
        SumContainer{MomentumSum{'q'}, IndexSum{ Index::Sigma }},
        std::vector<Operator>({
            c_q_sigma.hermitian_conjugate(),
            c_q_sigma
        })
    );

    const Term H_int(1, 
        Coefficient("U", MomentumList({Momentum('q'), Momentum('p'), Momentum('r')})),
        SumContainer{MomentumSum{'q', 'p', 'r'}, IndexSum{}},
        std::vector<Operator>({
            c_q_up.hermitian_conjugate(),
            c_q_down.with_momentum('p').hermitian_conjugate(),
            c_q_down.with_momentum(Momentum("p-r")),
            c_q_up.with_momentum(Momentum("q+r"))
        })
    );

    const std::vector<Term> hamiltonian = {H_kin, H_int};

    const std::vector<WickOperatorTemplate> wick_templates({ 
        WickOperatorTemplate({Num_Comparison}, Momentum(), OperatorType::Number) 
    });

    experimental::WickOrderedCollector hamiltonian_result = experimental::wick_decompose(hamiltonian, wick_templates);
    experimental::clean_wick_ordered_terms(hamiltonian_result);
    
    std::cout << "\\begin{align*}\n\t" << H_kin << "\\phantom{=}&\\\\\n" << H_int << "="
        << hamiltonian_result 
        << "\\end{align*}" << std::endl;

    const Term hexatic(1, std::vector<Operator>({
        c_q_up.hermitian_conjugate(),
        c_q_down.with_momentum('p').hermitian_conjugate(),
        c_q_down.with_momentum('k').hermitian_conjugate(),
        c_q_down.with_momentum('k'),
        c_q_down.with_momentum(Momentum("p")),
        c_q_up.with_momentum(Momentum("q"))
    }));
    experimental::WickOrderedCollector hexatic_result = experimental::wick_decompose(hexatic, wick_templates);
    experimental::clean_wick_ordered_terms(hexatic_result);

    std::cout << "\\begin{align*}\n\t" << hexatic << "=" << hexatic_result << "\\end{align*}" << std::endl;

    return 0;
}
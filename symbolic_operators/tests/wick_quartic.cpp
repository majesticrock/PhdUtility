#include <mrock/symbolic_operators/ExpectationValues>

#include <vector>

using namespace mrock::symbolic_operators;

int main() {
    const std::vector<WickOperatorTemplate> templates(
        {WickOperatorTemplate{{Num_Comparison}, Momentum(), OperatorType::Number},
         WickOperatorTemplate{{SC_Comparison}, Momentum(), OperatorType::SC}});

    const Momentum k('k');
    const Momentum p('p');
    const Operator c_k{k, Index::SpinUp, false};
    const Operator c_p{p, Index::SpinUp, false};
    const Operator c_k_dagger{k, Index::SpinUp, true};
    const Operator c_p_dagger{p, Index::SpinUp, true};
    const Operator c_minus_k{-k, Index::SpinDown, false};
    const Operator c_minus_p{-p, Index::SpinDown, false};
    const Operator c_minus_p_dagger{-p, Index::SpinDown, true};
    const Operator c_minus_k_dagger{-k, Index::SpinDown, true};

    const TermCollector density_density({Term(1, std::vector<Operator>{c_k_dagger, c_p_dagger, c_p, c_k})});
    const TermCollector pair_channel({Term(1, std::vector<Operator>{c_k_dagger, c_minus_k_dagger, c_minus_p, c_p})});
    const TermCollector annihilation_then_creation(
        {Term(1, std::vector<Operator>{c_minus_k, c_k, c_minus_p_dagger, c_p_dagger})});
    const TermCollector interleaved({Term(1, std::vector<Operator>{c_k_dagger, c_p, c_p_dagger, c_k})});

    const WickTermCollector density_wicks = wicks_theorem(density_density, templates, true);
    const WickTermCollector pair_wicks = wicks_theorem(pair_channel, templates, true);
    const WickTermCollector annihilation_then_creation_wicks =
        wicks_theorem(annihilation_then_creation, templates, true);
    const WickTermCollector interleaved_wicks = wicks_theorem(interleaved, templates, true);

    std::cout << "% Wick's theorem for a density-density quartic term\n\\begin{align*}\n\t" << "\\langle"
              << density_density << "\\rangle = " << density_wicks << "\\end{align*}" << std::endl;
    std::cout << "% Wick's theorem for a pair-channel quartic term\n\\begin{align*}\n\t" << "\\langle" << pair_channel
              << "\\rangle = " << pair_wicks << "\\end{align*}" << std::endl;
    std::cout << "% Wick's theorem for $cc c^{\\dagger} c^{\\dagger}$\n\\begin{align*}\n\t" << "\\langle"
              << annihilation_then_creation << "\\rangle = " << annihilation_then_creation_wicks << "\\end{align*}"
              << std::endl;
    std::cout << "% Wick's theorem for $c^{\\dagger} c c^{\\dagger} c$\n\\begin{align*}\n\t" << "\\langle"
              << interleaved << "\\rangle = " << interleaved_wicks << "\\end{align*}" << std::endl;

    WickTerm density_exchange;
    density_exchange.multiplicity = -1;
    density_exchange.delta_momenta.push_back({p, k});
    density_exchange.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);
    density_exchange.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);

    WickTerm density_direct;
    density_direct.multiplicity = 1;
    density_direct.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);
    density_direct.operators.emplace_back(OperatorType::Number, false, p, Index::SpinUp);

    WickTerm pair_condensate;
    pair_condensate.multiplicity = 1;
    pair_condensate.operators.emplace_back(OperatorType::SC, false, p);
    pair_condensate.operators.emplace_back(OperatorType::SC, true, k);

    WickTerm pair_density;
    pair_density.multiplicity = 1;
    pair_density.delta_momenta.push_back({p, k});
    pair_density.operators.emplace_back(OperatorType::Number, false, -k, Index::SpinDown);
    pair_density.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);

    WickTerm reversed_pair_holes;
    reversed_pair_holes.multiplicity = -1;
    reversed_pair_holes.delta_momenta.push_back({p, k});
    reversed_pair_holes.operators.emplace_back(OperatorType::HoleNumber, false, -k, Index::SpinDown);
    reversed_pair_holes.operators.emplace_back(OperatorType::HoleNumber, false, k, Index::SpinUp);

    WickTerm reversed_pair_sc;
    reversed_pair_sc.multiplicity = -1;
    reversed_pair_sc.operators.emplace_back(OperatorType::SC, false, k);
    reversed_pair_sc.operators.emplace_back(OperatorType::SC, true, p);

    WickTerm interleaved_density;
    interleaved_density.multiplicity = 1;
    interleaved_density.delta_momenta.push_back({p, k});
    interleaved_density.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);
    interleaved_density.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);

    WickTerm interleaved_hole;
    interleaved_hole.multiplicity = 1;
    interleaved_hole.operators.emplace_back(OperatorType::Number, false, k, Index::SpinUp);
    interleaved_hole.operators.emplace_back(OperatorType::HoleNumber, false, p, Index::SpinUp);

    WickTermCollector expected_density_wicks{density_exchange, density_direct};
    WickTermCollector expected_pair_wicks{pair_density, pair_condensate};
    const WickTermCollector expected_annihilation_then_creation_wicks{reversed_pair_holes, reversed_pair_sc};
    const WickTermCollector expected_interleaved_wicks{interleaved_density, interleaved_hole};
    if (density_wicks.terms != expected_density_wicks.terms || density_wicks[0].multiplicity != -1 ||
        density_wicks[1].multiplicity != 1 || pair_wicks.terms != expected_pair_wicks.terms ||
        pair_wicks[0].multiplicity != 1 || pair_wicks[1].multiplicity != 1 ||
        annihilation_then_creation_wicks.terms != expected_annihilation_then_creation_wicks.terms ||
        interleaved_wicks.terms != expected_interleaved_wicks.terms) {
        return 1;
    }
    return 0;
}
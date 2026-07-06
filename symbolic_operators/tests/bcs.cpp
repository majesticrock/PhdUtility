#include <cmath>
#include <numbers>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <map>
#include <stdexcept>
#include <iomanip>

#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/Wick.hpp>

using namespace mrock::symbolic_operators;

constexpr int N = 1000;

/* We consider a simple 1D chain with nearest-neighbor hopping => epsilon(k) = 2 cos(k) */
double epsilon(double k) {
    return -2. * std::cos(k);
}

/* 
We consider a BCS interaction g(q,p) = g Theta(omega - |epsilon(q)|) Theta(omega - |epsilon(p)|)
It therefore actually only depends on epsilon(q) and epsilon(p), which we leverage
We set g=0.1 and omega=0.2 */
constexpr double interaction_strength = 0.1 / static_cast<double>(N);
constexpr double omega = 0.2;
double g(double epsilon_q, double epsilon_p) {
    if (std::abs(epsilon_q) >= omega) {
        return double{};
    }
    if (std::abs(epsilon_p) >= omega) {
        return double{};
    }
    return interaction_strength;
};

/* 
Define the (mean-field) expectation values.
They depend on the order parameter Delta and the single-particle energy epsilon */
double expec_number(double energy, double delta) {
    if (std::abs(energy) >= omega) {
        return static_cast<double>(std::signbit(energy));
    }
    return 0.5 * (1. - energy / std::sqrt(energy*energy + delta*delta));
};
double expec_pair(double energy, double delta) {
    if (std::abs(energy) >= omega) {
        return 0.0;
    }
    return -0.5 * delta / std::sqrt(energy*energy + delta*delta);
};

/* Evaluates the momentum expression */
double evaluate_momentum(const Momentum& expression, const std::map<MomentumSymbol::name_type, double>& momentum_map) {
    double value{};
    for (const MomentumSymbol& momentum : expression) {
        value += momentum.factor * momentum_map.at(momentum.name);
    }
    if (expression.add_Q) {
        /* 
        Q is a special symbol with the property that 2Q is a reciprocal lattice vector.
        => in this model: Q=pi */
        value += std::numbers::pi;
    }
    return value;
}

/* Evaluates the 2 types of coefficients, also works with strings of coefficients */
double evaluate_coefficients(const std::vector<Coefficient>& coefficients, const std::map<MomentumSymbol::name_type, double>& momentum_map) {
    double value{1.0};
    for (const auto& coeff : coefficients) {
        if (coeff.name == "\\epsilon") {
            double momentum_value  = evaluate_momentum(coeff.momenta[0], momentum_map);
            value *= epsilon(momentum_value);
        }
        else if (coeff.name == "g") {
            double first_momentum  = evaluate_momentum(coeff.momenta[0], momentum_map);
            double second_momentum = evaluate_momentum(coeff.momenta[1], momentum_map);
            value *= g(epsilon(first_momentum), epsilon(second_momentum));
        }
        else {
            std::runtime_error("Coefficient not recognized!");
        }
    }
    return value;
}

/* Evaluates a string of WickOperator objects */
double evaluate_wick_operators(const std::vector<WickOperator>& operators, const std::map<MomentumSymbol::name_type, double>& momentum_map, double delta) {
    double value{1.0};
    for (const auto& op : operators) {
        double momentum_value = evaluate_momentum(op.momentum, momentum_map);

        switch(op.type) {
            case OperatorType::SC:
                value *= expec_pair(epsilon(momentum_value), delta);
                break;
            case OperatorType::Number:
                value *= expec_number(epsilon(momentum_value), delta);
                break;
            default:
                throw std::runtime_error("WickOperator not recognized!");
        }
    }
    return value;
}

/* 
The logic for computing an expression already uses that the overall structre is known. 
You can always do that in practice, since you can always just compute the Wick expressions (with the library)
and have a look at them. 
If you would rather implement a general handler, feel free to do so though. */
double evaluate_expression(const WickTerm& term, const std::array<double, N>& ks, double delta, double k=0) {
    /* The result contains either a single sum over q oder a double sum over p and q */
    std::map<MomentumSymbol::name_type, double> momentum_map;
    momentum_map['k'] = k;

    double value{};
    if (term.sums.momenta.size() == 1U) {
        for (int q=0; q<N; ++q) {
            momentum_map[term.sums.momenta[0]] = ks[q];
            value += evaluate_coefficients(term.coefficients, momentum_map) * evaluate_wick_operators(term.operators, momentum_map, delta);
        }
    }
    else if (term.sums.momenta.size() == 2U) {
        for (int q=0; q<N; ++q) {
            momentum_map[term.sums.momenta[0]] = ks[q];
            double psum{};
            for (int p=0; p<N; ++p) {
                momentum_map[term.sums.momenta[1]] = ks[p];
                psum += evaluate_coefficients(term.coefficients, momentum_map) * evaluate_wick_operators(term.operators, momentum_map, delta);
            }
            value += psum;
        }
    }
    else if (term.sums.momenta.empty()) {
        value = evaluate_coefficients(term.coefficients, momentum_map) * evaluate_wick_operators(term.operators, momentum_map, delta);
    }
    else {
        throw std::runtime_error("Number of momentum sums not supported!");
    }
    value *= static_cast<double>(term.multiplicity);
    // Account for spin degeneracy
    value *= std::pow(2.0, term.sums.spins.size());
    return value;
};

int main(int argc, char** argv) {
    // Setup a Momentum object representing 'k'
    const Momentum base_k = Momentum('k');

    // Setup a few operators to be used later
    const Operator c_k = Operator{ base_k, Index::SpinUp, false };
    const Operator c_minus_k = Operator{ -base_k, Index::SpinDown, false };
    const Operator c_k_dagger = Operator{ base_k, Index::SpinUp, true };
    const Operator c_minus_k_dagger = Operator{ -base_k, Index::SpinDown, true };

    // Setup a simple BCS Hamiltonian
    // Kinetic part: sum_(q,sigma) epsilon (q) c_(k,sigma)^dagger c_(k,sigma)
    const Term H_Kin(1, /* The prefactor of the term is 1 */
        Coefficient("\\epsilon", Momentum('q')), /* The coefficient of the term is epsilon (q) */
        SumContainer{ MomentumSum({ 'q' }), Index::Sigma }, /* The term contains a sum over the momentum q and the index sigma */
	    std::vector<Operator>({ /* The term has 2 operators: c_(q,sigma)^dagger and c_(q,sigma) */
		    Operator('q', 1, false, Index::Sigma, true), Operator('q', 1, false, Index::Sigma, false)
		}));

    // Pairing part: sum_(q,p) g(q, p) c_(q,up)^dagger c_(-q,down)^dagger c_(-p,down) c_(p,up)
    // The coefficent g(q, p) has the following symmetries:
    // g(q, p) = g(-q, p) = g(-q, -p) = g(q, -p) = g(p, q) = g^* (q, p)
    const Term H_Ph(-1, /* The prefactor of the term is -1 */
        Coefficient::RealInversionSymmetric("g", /* The coefficent is named 'g', 
            it is real g=g^* and inversion symmetric g(p) = g(-p) */
            MomentumList({ 'q', 'p' }), // The coefficient depends on the momenta q and p
            std::function<void(Coefficient&)>([](Coefficient& coeff) {  /* The coefficient has an additional custom symmetry
                that is being handled by this lamdba expression. In this case, g(p, q) = g(q, p)
                so the lambda just sorts the momenta to bring all coefficients to the same notation. */
                coeff.momenta.sort(); 
            })
        ),
		SumContainer{ MomentumSum({ 'p', 'q' }) },  /* The term contains a sum over the momenta q and p*/
		std::vector<Operator>({ /* The term has 4 operators: c_(q,up)^dagger c_(-q,down)^dagger c_(-p,down) c_(p,up) */
			c_k_dagger.with_momentum('q'), c_minus_k_dagger.with_momentum('q'),
			c_minus_k.with_momentum('p'), c_k.with_momentum('p') })
		);

    // The Hamiltonian is represented by a vector of the two Term objects above
    const std::vector<Term> H{ H_Kin, H_Ph };

    // Lets say, we have a superconducting system, so for applying Wick's theorem, 
    // we focus on number and pair creation/annihilation operators
    const std::vector<WickOperatorTemplate> templates({
            /* Template for pair creation/annhilation operators */
			WickOperatorTemplate{ 
                {IndexComparison{ /* Alternatively, the constexpr object SC_Comparison could be used */
                    false, /* The indizes (here only spins), must be fixed */
                    Index::SpinUp, /* The second index must be SpinUp */
                    Index::SpinDown} /* The first index must be SpinDown */
                },  
                /* Note that the pair annihilation operator is defined as f_k := c_(-k,down) c_(k,up).
                Therefore, the 'base operatore' is the second one in the expression c_(k,up),
                and we have to give the index of the second operator first. */
                Momentum(), /* The total momentum of the expression must be 0 */
                OperatorType::SC /* It is an OperatorType::SC term, i.e., cc or c^dagger c^dagger with total Momentum 0 */
            },
            /* Template for number operators */
			WickOperatorTemplate{ 
                {IndexComparison{ /* Alternatively, the constexpr object Num_Comparison could be used */
                    true }  /* The indizes (here: spins) of both operators must be equal, but no other restriction is placed */
                },
                Momentum(), /* The total momentum of the expression must be 0 */
                OperatorType::Number /* It is an OperatorType::Number term, i.e., c^dagger c with total Momentum 0 */
            }
		});
    
    /* Name the symmetries that the expectation values have */
    std::vector<std::unique_ptr<WickSymmetry>> symmetries;
    /* The expectation values do not depends on the spin, i.e., f(down) = f(up) */
    symmetries.push_back(std::make_unique<SpinSymmetry>());
    /* The expectation values are inversion symmetric, i.e., f(k) = f(-k) */
    symmetries.push_back(std::make_unique<InversionSymmetry>());
    /* The expectation values of the SC-type are real, i.e., f = f^* */
    symmetries.push_back(std::make_unique<PhaseSymmetry<OperatorType::SC>>());

    /* 
    Let us do something quantitative.
    We want to compute the expectation of the Hamiltonian <H>.
    To do so, we must define g(q,p), epsilon(q), and the expecation values 
    of the number operators <n_k> and of the pair creation/annihilation operators <f_k> */

    
    /* We precompute all epsilon and for all k*/
    constexpr double dk = 2. * std::numbers::pi / static_cast<double>(N);
    std::array<double, N> ks{};
    std::generate(ks.begin(), ks.end(), [i=0]() mutable {
        return (i++ * dk) - std::numbers::pi;
    });
    std::array<double, N> epsilons{};
    std::transform(ks.begin(), ks.end(), epsilons.begin(), [](double k) {
        return epsilon(k);
    });

    /* We precompute all  */
    /* To find Delta, we solve the mean-field self-consistency problem */
    auto self_consistency = [&epsilons](double delta) {
        double new_delta{};
        for (auto eps : epsilons) {
            if (std::abs(eps) < omega) {
                new_delta += 1. / std::sqrt(eps*eps + delta*delta);
            }
        }
        new_delta *= 0.5 * interaction_strength * delta;
        return new_delta;
    };
    double delta = 0.1;
    double new_delta = 0.1;
    do {
        delta = new_delta;
        new_delta = self_consistency(delta);
    } while (std::abs(delta - new_delta) > 1e-4);
    std::cout << "Self-consistency problem solved with Delta=" << new_delta << std::endl;

    /* Compute the expecation value of H with the analytical formula for later comparison */
    /* Kinetic part, the factor of 2 accounts for the spins */
    const double kinetic = 2 * std::accumulate(epsilons.begin(), epsilons.end(), double{}, 
        [&delta](const double& current, const double& epsilon) {
            return current + epsilon * expec_number(epsilon, delta);
        }
    );
    /* Pairing part */
    const double pairing = std::accumulate(epsilons.begin(), epsilons.end(), double{},
        [&delta](const double& current, const double& epsilon) {
            if (std::abs(epsilon) >= omega) {
                return current;
            }
            return current + expec_pair(epsilon, delta);
        }
    );
    /* The <n_k><n_k> contraction; it scales as 1/N, and can be neglected in the thermodynamic limit. */
    const double one_over_n_part = std::accumulate(epsilons.begin(), epsilons.end(), double{},
        [&delta](const double& current, const double& epsilon) {
            if (std::abs(epsilon) >= omega) {
                return current;
            }
            return current + expec_number(epsilon, delta) * expec_number(epsilon, delta);
        }
    );
    const double analytical = (kinetic - interaction_strength * (pairing*pairing + one_over_n_part)) / N;
    /* 
    In addition to those terms, there is a <n><n> contraction of the interaction term
    However, this one scales as 1/N and may thus be neglected in the thermodynamic limit. */
    
    /* Now compute the expecation value using the implementation of Wick's theorem */
    WickTermCollector wicks;
    wicks_theorem(H, templates, wicks);
    /* Clean up the result and apply symmetries */
    clean_wicks(wicks, symmetries);
    std::cout << "Expressional result from Wick's theorem:\n" << wicks << "\n" << std::endl;

    double numerical{};
    for (const auto& term : wicks) {
        numerical += evaluate_expression(term, ks, delta);
    }
    numerical /= N;
    std::cout << std::setprecision(14) 
              << "Analytical value for <H>/N = " << analytical << std::endl;
    std::cout << " Numerical value for <H>/N = " << numerical  << std::endl;
    if (std::abs(numerical - analytical) > 1e-10) {
        return -1;
    }

    ///////////////////////////////////////////////////////////////
    /* Next we evaluate <[H, f_k]>, where f_k is a pair creation operator */
    Term right(1, std::vector<Operator>({ c_minus_k, c_k }));

    /* Compute the commutator and clean up the result */
    std::vector<Term> commutator_result = commutator(H, right);
    clean_up(commutator_result);
    std::cout << "Result of the commutator:\n" << commutator_result << "\n" << std::endl;

    /* Apply Wick's theorem and clean up the result */
    WickTermCollector commutator_wicks;
    wicks_theorem(commutator_result, templates, commutator_wicks);
    clean_wicks(commutator_wicks, symmetries);
    std::cout << "Expressional result from Wick's theorem:\n" << commutator_wicks << "\n" << std::endl;


    /* Compute the expecation value of <[H,f_k]> with the analytical formula for later comparison */
    /* Kinetic part, the factor of 2 accounts for the spins */

    const int k = N/4 - 2;
    double bilinear = std::accumulate(epsilons.begin(), epsilons.end(), double{},
        [&delta](const double& current, const double& epsilon) {
            if (std::abs(epsilon) >= omega) {
                return current;
            }
            return current + expec_pair(epsilon, delta);
        }
    );
    bilinear *= -interaction_strength;
    if (std::abs(epsilons[k]) < omega) {
        bilinear += 2. * epsilons[k] * expec_pair(epsilons[k], delta);
    }

    /* Pairing part */
    double quartic = std::accumulate(epsilons.begin(), epsilons.end(), double{},
        [&delta](const double& current, const double& epsilon) {
            if (std::abs(epsilon) >= omega) {
                return current;
            }
            return current + expec_pair(epsilon, delta);
        }
    );
    quartic *= 2. * expec_number(epsilons[k], delta) * interaction_strength;

    if (std::abs(epsilons[k]) >= omega) {
        quartic = 0.0;
        bilinear = 0.0;
    }
    const double analytical_commutator = (bilinear + quartic) ;// N;

    double numerical_commutator{};
    for (const auto& term : commutator_wicks) {
        numerical_commutator += evaluate_expression(term, ks, delta, ks[k]);
    }
    //numerical_commutator /= N;
    std::cout << std::setprecision(14) 
              << "Analytical value for <[H, f_k]>/N = " << analytical_commutator << std::endl;
    std::cout << " Numerical value for <[H, f_k]>/N = " << numerical_commutator  << std::endl;
    if (std::abs(numerical_commutator - analytical_commutator) > 1e-10) {
        return -1;
    }


    return 0;
}
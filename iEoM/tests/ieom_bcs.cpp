#include <mrock/iEoM/XPResolvent.hpp>
#include <cmath>
#include <numbers>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <array>

using namespace mrock::iEoM;

/**
 * Assumed operator basis:
 * 
 * S_k^+ =  {   
 *              f_k + f_k^+dagger
 *              n_(k,up) + n_(-k down)
 *          }
 * S_k^- =  {   
 *              f_k - f_k^+dagger
 *          }
 * 
 * inversion symmetry is assumed
 */

constexpr int N_k = 4;
constexpr int N_lanczos = N_k / 10;
constexpr double U = -0.5;
std::array<double, N_k> epsilons;

struct BCSTester : public XPResolvent<double> {
    double Delta;
    
    double E(double eps) {
        return std::sqrt(eps*eps + Delta*Delta);
    }
    double number(double eps) {
        return 0.5 * (1. - eps / E(eps));
    }
    double pair(double eps) {
        return -0.5 * Delta / E(eps);
    }
    
    void fill_M() {
        K_plus.resize(2*N_k, 2*N_k);
        K_minus.resize(N_k, N_k);

        for (int k=0; k<N_k; ++k) {
            // K-diagonal filling
            double number_k = number(epsilons[k]);
            double pair_k = pair(epsilons[k]);

            // The terms with epsilon actually read epsilon - mu + U/2
            // But at half-filling mu = U/2
            K_plus(k*2, k*2) = 4 * epsilons[k] * (1 - 2*number_k);
            K_plus(k*2 + 1, k*2 + 1) = -8 * Delta * pair_k;

            K_plus(k*2 + 1, k*2) = -8 * epsilons[k] * pair_k;
            K_plus(k*2, k*2 + 1) = 4 * Delta * ( 1 - 2*number_k );

            K_minus(k, k) = 4 * epsilons[k] * (1 - 2*number_k ) - 8 * Delta * pair_k;

            for (int l=0; l<N_k; ++l) {
                // K-offdiagonal filling, these terms tend to be small in
                // the thermodynamic limit, as they scale as 1/N - but there
                // are a lot of them. In fact, these terms are what distinguishes
                // the method from standard mean-field. They carry collective behavior
                double number_l = number(epsilons[l]);
                double pair_l = pair(epsilons[l]);

                K_plus(k*2, l*2) += (2 * U / N_k) * (1 - 2*number_k) * (1 - 2*number_l);
                K_plus(k*2 + 1, l*2 + 1) += (8 * U / N_k) * pair_l * pair_k;

                K_plus(k*2 + 1, l*2) += (-4 * U / N_k) * pair_l * (1 - 2*number_k);
                K_plus(k*2, l*2 + 1) += (-4 * U / N_k) * pair_k * (1 - 2*number_l);

                K_minus(k, l) += (2 * U / N_k) * (1 - 2*number_k) * (1 - 2*number_l) + (8 * U / N_k) * pair_l * pair_k;
            }
        }
    }

    void fill_matrices() {
        fill_M();
        L.resize(2 * N_k, N_k);

        for (int k=0; k<N_k; ++k) {
            // L has only k-diagonal components
            L(k, k) = -2 * (1 - 2*number(epsilons[k]));
            L(k+1, k) = 4 * pair(epsilons[k]);
        }
    }

    void create_starting_states() {
        starting_states.clear();
        // Vector is defined as a protected member of the algorithm.
        // Here, it is an Eigen::VectorXd.
        starting_states.push_back(StartingState<double>{  Vector::Ones(N_k), Vector::Ones(2*N_k), "SC" });
    }

    BCSTester(double _delta) 
        : XPResolvent<double>(1e-6, 2*N_k, N_k),
        Delta{_delta} {}
};

void save_data(const std::vector<BCSTester::ResolventReturnData>& data, std::ofstream& out)
{
	out << std::scientific << std::setprecision(10);
    std::array<std::string, N_lanczos+1> lines;
    lines[0] = "#";

	for (const auto& data_wrapper : data) {
        lines[0] += data_wrapper.name + ":a_i " + data_wrapper.name + ":b_i ";

		for (const auto& resolvent_data : data_wrapper.lanczos) {
            for (int i=0; i<N_lanczos; ++i) {
                lines[i+1] += " " + std::to_string(resolvent_data.a_i[i]) + " " + std::to_string(resolvent_data.b_i[i]);
            }
        }
	}

    for (const auto& line : lines) {
        out << line << '\n';
    }
};

int main() {
    for (int k=0; k<N_k; ++k) {
        epsilons[k] += -2. * std::cos(std::numbers::pi * (1. - static_cast<double>(k) / static_cast<double>(N_k)) );
    }

    auto selfconsistency = [](double delta_old) {
        double delta_new{};
        for (int k=0; k<N_k; ++k) {
            delta_new += -0.5 * delta_old / std::sqrt(epsilons[k]*epsilons[k] + delta_old*delta_old);
        }
        delta_new *= U/N_k;
        return delta_new;
    };

    double delta_old = 0.1;
    double delta_new = selfconsistency(delta_old);
    while(std::abs(delta_old - delta_new) > 1e-6) {
        delta_old = delta_new;
        delta_new = selfconsistency(delta_old);
    }
    std::cout << std::scientific << "Self-consistency converged! Delta = " << delta_new << std::endl; 

    BCSTester tester(delta_new);
    std::vector<BCSTester::ResolventReturnData> result = tester.compute_collective_modes<10>(N_lanczos);

    std::ofstream filestream("bcs_result.txt", std::ios_base::out);
    if (filestream.is_open()) {
        save_data(result, filestream);
        filestream.close();
    }
    else {
        std::cerr << "ieom_bcs.cpp could not open output filestream" << std::endl;
        return 1;
    }
    std::cout << "Test data has been saved to bcs_result.txt" << std::endl;

    return 0;
}
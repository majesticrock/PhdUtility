#include <fstream>
#include <sstream>
#include <Utility/LaTeXOutput.hpp>
#include "Definitions/Hubbard.hpp"
#include "Definitions/HubbardDispersions.hpp"
#include "Definitions/Continuum.hpp"
#include <SymbolicOperators/Wick.hpp>
#include <memory>
#include <filesystem>

using namespace SymbolicOperators;
using term_vec = std::vector<Term>;
using op_vec = std::vector<Operator>;

std::unique_ptr<StandardOperators> get_model(std::string const& model_type) {
	if (model_type == "hubbard") {
		return std::make_unique<Hubbard>();
	}
	else if (model_type == "continuum") {
		return std::make_unique<Continuum>();
	}
	else if (model_type == "hubbard_dispersions") {
		return std::make_unique<HubbardDispersions>();
	}
	else {
		throw std::invalid_argument("Model not recognized! " + model_type);
	}
}

int main(int argc, char** argv) {
	const std::string save_folder = "../commutators/";
	//WickTerm parse_test("1 sum:momentum{p,q} c:V{p;q} o:n{k-p-3x;up} o:f{k+l;}");
	//std::cout << parse_test << std::endl;

	constexpr bool print = true;
	if (argc < 3) {
		std::cerr << "Syntax: ./build/main <XP/std> <model>" << std::endl;
		return 1;
	}
	const std::string EXECUTION_TYPE = argv[1];
	const std::string MODEL_TYPE = argv[2];

	if (EXECUTION_TYPE == "test") {
		WickTerm wick;
		wick.multiplicity = 1;
		wick.temporary_operators = { Hubbard::c_minus_k_Q, Hubbard::c_k_Q,
			Hubbard::c_k_Q_dagger, StandardOperators::c_k };
		auto wick_results = identifyWickOperators(wick, Hubbard().templates());
		std::cout << "Testing on: $" << wick.temporary_operators << "$\n\n";
		std::cout << "Pre clean:\n\n" << Utility::as_LaTeX(wick_results, "align*") << std::endl;
		cleanWicks(wick_results, Hubbard().symmetries());
		std::cout << "Post clean:\n\n" << Utility::as_LaTeX(wick_results, "align*") << std::endl;

		return 0;
	}

	const std::string name_prefix = EXECUTION_TYPE == "XP" ? "XP_" : "";
	const bool debug = EXECUTION_TYPE == "debug";

	const std::unique_ptr<StandardOperators> model = get_model(MODEL_TYPE);
	const std::string sub_folder = model->get_subfolder();
	if (!debug)
		std::filesystem::create_directories(save_folder + sub_folder);

	const term_vec H = model->hamiltonian();
	const std::vector<WickOperatorTemplate> templates = model->templates();

	std::vector<term_vec> basis;
	if (EXECUTION_TYPE == "XP") {
		basis = model->XP_basis();
	}
	else if (EXECUTION_TYPE == "std") {
		basis = model->STD_basis();
	}
	else if (debug) {
		basis = {
			// 0: f + f^+
			std::vector<Term>({
				Term(1, std::vector<Operator>({ Operator(Momentum({{-1, 'k'}, {-1, 'x'}}), SpinDown, false), StandardOperators::c_k}))
			})
		};
	}
	const auto symmetries = model->symmetries();

	std::vector<term_vec> basis_daggered(basis);
	for (auto& t : basis_daggered) {
		hermitianConjugate(t);
		renameMomenta(t, 'k', 'l');
		if (debug) {
			renameMomenta(t, 'x', 'y');
		}
	}
	std::cout << "\\begin{align*}\n\t H =" << H << "\\end{align*}" << std::endl;
	for (size_t i = 0U; i < basis.size(); ++i)
	{
		term_vec commute_with_H;
		commutator(commute_with_H, H, basis[i]);
		cleanUp(commute_with_H);
		if (debug)
			std::cout << "\\begin{align*}\n\t[ H, " << toStringWithoutPrefactor(basis[i]) << " ] ="
			<< commute_with_H << "\\end{align*}" << std::endl;

		for (size_t j = 0U; j < basis.size(); ++j)
		{
			if (print) {
				std::cout << "\\subsection{" << i << "." << j << "}" << std::endl;
			}
			term_vec terms;
			commutator(terms, basis_daggered[j], commute_with_H);
			cleanUp(terms);

			if (debug)
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
				<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << terms << "\\end{align*}" << std::endl;

			WickTermCollector wicks;
			wicks_theorem(terms, templates, wicks);
			clearEtas(wicks);
			cleanWicks(wicks, symmetries);

			for(auto& wickterm : wicks){
				if(wickterm.coefficients.front().name == "\\rho"){
					wickterm.sums.push_back(Index::SigmaPrime);
					wickterm.sums.push_back('r');
					wickterm.coefficients.front() = Coefficient::parse_string("V{0;}");
					wickterm.operators.push_back(WickOperator("n{r;sigma'}"));
				}
			}
			cleanWicks(wicks, symmetries);

			if (debug || print) {
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
					<< ", [H, " << toStringWithoutPrefactor(basis[i]) << " ]] =" << wicks << "\\end{align*}" << std::endl;
			}

			// serialization
			if (!debug) {
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs(save_folder + sub_folder + name_prefix + "wick_M_" + std::to_string(j) + "_" + std::to_string(i) + ".txt");
				boost::archive::text_oarchive oa(ofs);
				oa << wicks;
				ofs.close();
			}

			terms.clear();
			wicks.clear();
			commutator(terms, basis_daggered[j], basis[i]);
			cleanUp(terms);
			wicks_theorem(terms, templates, wicks);
			clearEtas(wicks);
			cleanWicks(wicks, symmetries);

			if (debug || print)
				std::cout << "\\begin{align*}\n\t[ " << toStringWithoutPrefactor(basis_daggered[j])
				<< ", " << toStringWithoutPrefactor(basis[i]) << " ] =" << wicks << "\\end{align*}" << std::endl;
			// serialization
			if (!debug) {
				// create an output file stream and a text archive to serialize the vector
				std::ofstream ofs(save_folder + sub_folder + name_prefix + "wick_N_" + std::to_string(j) + "_" + std::to_string(i) + ".txt");
				boost::archive::text_oarchive oa(ofs);
				oa << wicks;
				ofs.close();
			}
		}
	}

	return 0;
}

/* Output code if needed
*
Term right(1, Coefficient(), op_vec({
	c_minus_k, c_k
	}));
Term left(1, Coefficient(), op_vec({
	c_l_dagger, c_minus_l_dagger
	}));

std::cout << "\\begin{align*}\n\t\\langle [" << left.toStringWithoutPrefactor() << ", " << right.toStringWithoutPrefactor() << "] \\rangle = " << wicks << "\\end{align*}" << std::endl;
std::cout << "\\begin{align*}\n\t\\langle [" << left.toStringWithoutPrefactor() << ", [ H, " << right.toStringWithoutPrefactor() << "] ] \\rangle = " << wicks << "\\end{align*}" << std::endl;
*/

/* Example of how to to read the output
// create an input file stream and a text archive to deserialize the vector
std::ifstream ifs("wick_terms.txt");
boost::archive::text_iarchive ia(ifs);
WickTermCollector deserialized_terms;
ia >> deserialized_terms;
ifs.close();
*/
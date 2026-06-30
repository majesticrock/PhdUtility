#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/KroneckerDeltaUtility.hpp>
#include <mrock/utility/RangeUtility.hpp>
#include <sstream>

namespace mrock::symbolic_operators {
	Term::Term(IntFractional _multiplicity, std::vector<Coefficient> _coefficients, const SumContainer& _sums, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, _coefficients, _sums, _operators)
	{}

	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const SumContainer& _sums, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, _coefficient, _sums, _operators)
	{}

	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const MomentumSum& _sum_momenta, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, _coefficient, SumContainer{ _sum_momenta, {} }, _operators)
	{}

	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const IndexSum& _sum_indizes, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, _coefficient, SumContainer{ {}, _sum_indizes }, _operators)
	{}

	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, _coefficient, _operators)
	{}

	Term::Term(IntFractional _multiplicity, const SumContainer& _sums, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, std::vector<Coefficient>(), _sums, _operators)
	{}

	Term::Term(IntFractional _multiplicity, const MomentumSum& _sum_momenta, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, std::vector<Coefficient>(), SumContainer{ _sum_momenta, {} }, _operators)
	{}

	Term::Term(IntFractional _multiplicity, const IndexSum& _sum_indizes, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, std::vector<Coefficient>(), SumContainer{ {}, _sum_indizes }, _operators)
	{}

	Term::Term(IntFractional _multiplicity, const std::vector<Operator>& _operators)
		: AbstractTerm<Operator>(_multiplicity, _operators)
	{}

	void Term::print() const {
		std::cout << *this << std::endl;
	}

	bool Term::resolve_deltas()
	{
		if (!resolve_momentum_deltas()) return false;
		if (!resolve_index_deltas()) return false;

		// Check for the Pauli principle
		for (auto it = operators.begin(); it != operators.end(); ++it) {
			for (auto jt = it+1; jt != operators.end(); ++jt) {
				if (*it == *jt) return false; 
				if (it->is_daggered != jt->is_daggered) break;
			}
		}
		
		return true;
	}

	void Term::discard_zero_momenta() {
		for (auto& op : operators) {
			op.momentum.remove_zeros();
		}
		for (auto& coeff : coefficients) {
			coeff.momenta.remove_zeros();
		}
	}

	void Term::sort()
	{
		for (auto& coeff : coefficients) {
			for (auto& momentum : coeff.momenta) {
				momentum.sort();

				if (coeff.inversion_symmetry && !momentum.momentum_list.empty()) {
					if (momentum.momentum_list[0].factor < 0) {
						momentum.flip_momentum();
					}
				}
				if (coeff.Q_changes_sign && momentum.add_Q) {
					momentum.add_Q = false;
					flip_sign();
				}
			}
		}
		size_t new_n;
		size_t n = operators.size();
		while (n > 1U) {
			new_n = 0U;
			for (size_t i = 1U; i < n; ++i)
			{
				if (operators[i].is_daggered != operators[i - 1].is_daggered) continue;
				if (operators[i].is_fermion != operators[i - 1].is_fermion) continue;
				const Index l_idx = operators[i - 1].first_index();
				const Index r_idx = operators[i].first_index();
				if (operators[i].is_daggered) {
					// c^+ c^+
					if (r_idx == Index::SpinUp && l_idx != Index::SpinUp) {
						perform_operator_swap(operators[i], operators[i - 1]);
						new_n = i;
					}
					else if (l_idx == Index::SpinDown && r_idx != Index::SpinDown) {
						perform_operator_swap(operators[i], operators[i - 1]);
						new_n = i;
					}
					else if (l_idx > Index::SpinDown && r_idx > Index::SpinDown && l_idx > r_idx) {
						perform_operator_swap(operators[i], operators[i - 1]);
						new_n = i;
					}
				}
				else {
					// c c
					if (r_idx == Index::SpinDown && l_idx != Index::SpinDown) {
						perform_operator_swap(operators[i], operators[i - 1]);
						new_n = i;
					}
					else if (l_idx == Index::SpinUp && r_idx != Index::SpinUp) {
						perform_operator_swap(operators[i], operators[i - 1]);
						new_n = i;
					}
					else if(l_idx > Index::SpinDown && r_idx > Index::SpinDown && l_idx < r_idx) {
						perform_operator_swap(operators[i], operators[i - 1]);
						new_n = i;
					}
				}
			}
			n = new_n;
		}

		n = operators.size();
		while (n > 1U) {
			new_n = 0U;
			for (size_t i = 1U; i < n; ++i)
			{
				if (operators[i].is_daggered != operators[i - 1].is_daggered) continue;
				if (operators[i].first_index() != operators[i - 1].first_index()) continue;

				if(momentum_order(operators[i - 1].momentum, operators[i].momentum)) {
					perform_operator_swap(operators[i], operators[i - 1]);
					new_n = i;
				}
			}
			n = new_n;
		}

		// Sort the occurring coefficients in alphabetical order
		std::sort(coefficients.begin(), coefficients.end(), [](const Coefficient& a, const Coefficient& b) {
			if (a.name == b.name) {
				return (a.is_daggered && (!b.is_daggered));
			}
			return a.name < b.name;
		});

		// check whether we can swap the sign of each momentum in the coefficients
		// 26.04.2024, I have no idea what I did here, nor do I know why I did what I did
		for (const auto& coeff : coefficients) {
			if (!(coeff.inversion_symmetry)) return;
			if (std::any_of(coeff.momenta.begin(), coeff.momenta.end(), [](Momentum const& momentum) {
				return momentum.momentum_list.size() > 1U;
				})) return;
		}

		for (const auto& sum_mom : sums.momenta) {
			bool first_occurance = true;
			for (auto& op : operators) {
				int i = op.momentum.is_used_at(sum_mom);
				if (i > -1) {
					if (first_occurance) {
						if (op.momentum.momentum_list[i].factor < 0) {
							first_occurance = false;
						}
						else {
							break;
						}
					}
					op.momentum.momentum_list[i].factor *= -1;
				}
			}
		}
	}

	bool Term::is_equal(const Term& other) const {
		if (this->coefficients != other.coefficients) return false;
		if (this->sums != other.sums) return false;
		if (this->delta_indizes != other.delta_indizes) return false;
		if (this->delta_momenta != other.delta_momenta) return false;
		if (this->operators != other.operators) return false;
		return true;
	}

	bool Term::is_normal_ordered() const {
		for	(size_t i = 1U; i < operators.size(); ++i) {
			if (operators[i - 1].is_fermion == operators[i].is_fermion) {
				if(!operators[i - 1].is_daggered && operators[i].is_daggered) {
					return false;
				}
			}
		}
		return true;
	}

	std::string Term::to_string_without_prefactor() const
	{
		std::ostringstream os;
		if (!this->sums.spins.empty()) {
			os << "\\sum_{ ";
			for (const auto& index : this->sums.spins) {
				os << index << " ";
			}
			os << "}";
		}
		if (!this->sums.momenta.empty()) {
			os << "\\sum_{ ";
			for (const auto& momentum : this->sums.momenta) {
				os << momentum << " ";
			}
			os << "}";
		}
		os << this->coefficients << " ";
		for (const auto& delta : delta_momenta) {
			os << delta;
		}
		for (const auto& delta : delta_indizes) {
			os << delta;
		}

		if (this->is_identity()) {
			os << " \\hat{1} ";
			return os.str();
		}
		for (const auto& op : this->operators) {
			os << op << " ";
		}
		return os.str();
	}

	Term& Term::hermitian_conjugate_inplace() {
		std::reverse(this->operators.begin(), this->operators.end());
		for (auto& op : this->operators) {
			op.hermitian_conjugate_inplace();
		}
		for (auto& coeff : this->coefficients) {
			coeff.hermitian_conjugate_inplace();
		}
		return *this;
	}

	Term Term::hermitian_conjugate() const {
		Term copy(*this);
		copy.hermitian_conjugate_inplace();
		return copy;
	}

	void Term::rename_indizes(Index what, Index to) {
		if (what == to) return;
		for (auto& index_sum : sums.spins) {
			if (index_sum == to) {
				throw std::invalid_argument("You are replacing an index sum with an index that already exists!");
			}
			if (index_sum == what) {
				index_sum = to;
			}
		}
		for (auto& coeff : coefficients) {
			for (auto& index : coeff.indizes) {
				if(index == what) {
					index = to;
				}
			}
		}
		for (auto& op : operators) {
			for (auto& index : op.indizes) {
				if(index == what) {
					index = to;
				}
			}
		}
	}

	void Term::rename_momenta(const MomentumSymbol::name_type what, const MomentumSymbol::name_type to) {
		if (what == to) return;
		for (auto& mom_sum : sums.momenta) {
			if (mom_sum == to) {
				throw std::invalid_argument("You are replacing a momentum sum with an index that already exists!");
			}
			if (mom_sum == what) {
				mom_sum = to;
			}
		}
		for (auto& coeff : coefficients) {
			for (auto& mom : coeff.momenta) {
				mom.replace_occurances(what, Momentum(to));
			}
		}
		for (auto& op : operators) {
			op.momentum.replace_occurances(what, Momentum(to));
		}
	}

	void Term::swap_momenta(const MomentumSymbol::name_type a, const MomentumSymbol::name_type b) {
		this->rename_momenta(a, '_');
		this->rename_momenta(b, a);
		this->rename_momenta('_', b);
	}

	void Term::transform_momentum_sum(const MomentumSymbol::name_type what, const Momentum to, const MomentumSymbol::name_type new_sum_index) {
		auto pos = std::find(sums.momenta.begin(), sums.momenta.end(), what);
		if (pos == sums.momenta.end()) {
			throw std::invalid_argument("You are trying to perform a sum transformation on a momentum that is not being summed over!");
		} 
		else {
			*pos = new_sum_index;
		}
		for (auto& coeff : coefficients) {
			for (auto& mom : coeff.momenta) {
				mom.replace_occurances(what, to);
			}
		}
		for (auto& op : operators) {
			op.momentum.replace_occurances(what, to);
		}
	}

	void normal_order(std::vector<Term>& terms) {
		for (int t = 0; t < terms.size();) {
		normal_order_outerLoop:
			if (t >= terms.size()) break;
			size_t n = terms[t].operators.size();
			size_t new_n{};
			// First sort so that the bosons are upfront
			while (n > 1U) {
				new_n = 0U;
				for (size_t i = 1U; i < terms[t].operators.size(); ++i)
				{
					if (terms[t].operators[i - 1].is_fermion && !terms[t].operators[i].is_fermion) {
						new_n = i;
						std::swap(terms[t].operators[i - 1], terms[t].operators[i]);
					}
				}
				n = new_n;
			}
			
			n = terms[t].operators.size();
			new_n = 0U;
			while (n > 1U) {
				new_n = 0U;
				for (size_t i = 1U; i < terms[t].operators.size(); ++i)
				{
					if (!terms[t].operators[i - 1].is_fermion && terms[t].operators[i].is_fermion) continue;
					if (!(terms[t].operators[i - 1].is_daggered) && (terms[t].operators[i].is_daggered)) {
						bool other_deltas = false;
						new_n = i;
						// Swap cc^+
						std::swap(terms[t].operators[i - 1], terms[t].operators[i]);

						// Add a new term where cc^+ is replaced by the appropriate delta
						Term new_term(terms[t]);
						// flip the signs if we have fermions
						if (terms[t].operators[i - 1].is_fermion && terms[t].operators[i].is_fermion) {
							terms[t].flip_sign();
						}			
						if (new_term.operators[i - 1].indizes.size() != new_term.operators[i].indizes.size()) {
							throw std::invalid_argument("Operators do not have the same index count.");
						}

						if( (new_term.operators[i - 1].first_index() == Index::SpinUp && new_term.operators[i].first_index() == Index::SpinDown)
							|| (new_term.operators[i - 1].first_index() == Index::SpinDown && new_term.operators[i].first_index() == Index::SpinUp) ) {
								continue;
						}
						else if(new_term.operators[i - 1].first_index() != new_term.operators[i].first_index()) {
							new_term.delta_indizes.push_back(
								make_delta(new_term.operators[i - 1].first_index(), new_term.operators[i].first_index())
							);
						}
						for (int c = 1; c < new_term.operators[i - 1].indizes.size(); c++)
						{
							// if the indizes are not the same we emplace a delta
							// otherwise no action is required
							if (new_term.operators[i - 1].indizes[c] != new_term.operators[i].indizes[c]) {
								other_deltas = true;
								new_term.delta_indizes.push_back(
									make_delta(new_term.operators[i - 1].indizes[c], new_term.operators[i].indizes[c])
								);
							}
						}
						if (new_term.operators[i - 1].momentum != new_term.operators[i].momentum) {
							other_deltas = true;
							new_term.delta_momenta.push_back(
								make_delta(new_term.operators[i - 1].momentum, new_term.operators[i].momentum)
							);
						}
						else {
							other_deltas = true;
						}

						new_term.operators.erase(new_term.operators.begin() + i - 1, new_term.operators.begin() + i + 1);

						//if (new_term.resolve_deltas()) {
						if (other_deltas) terms.push_back(new_term);
						//}
					}
					else if (terms[t].operators[i - 1] == terms[t].operators[i]) {
						if (terms[t].operators[i - 1].is_fermion) {
							// two identical fermion operators = 0
							terms.erase(terms.begin() + t);
							goto normal_order_outerLoop;
						}
					}
				}
				n = new_n;
			}
			++t;
		}
	}

#define fill_reciever(x) reciever[0].x = left.x; mrock::utility::append_vector(reciever[0].x, right.x); reciever[1].x = left.x; mrock::utility::append_vector(reciever[1].x, right.x);
	std::vector<Term> commutator(const Term& left, const Term& right)
	{
		std::vector<Term> reciever(2);
		reciever[0] = left;
		reciever[0].multiplicity *= right.multiplicity;
		mrock::utility::append_vector(reciever[0].operators, right.operators);
		reciever[1] = right;
		reciever[1].multiplicity *= left.multiplicity;
		mrock::utility::append_vector(reciever[1].operators, left.operators);
		reciever[1].flip_sign();

		fill_reciever(coefficients);
		fill_reciever(sums.momenta);
		fill_reciever(sums.spins);
		fill_reciever(delta_momenta);
		fill_reciever(delta_indizes);

		normal_order(reciever);
		return reciever;
	}

	std::vector<Term> commutator(const std::vector<Term>& left, const std::vector<Term>& right)
	{
		std::vector<Term> reciever;
		reciever.reserve(2 * left.size() * right.size());
		
		for (const auto& left_term : left)
		{
			for (const auto& right_term : right)
			{
				std::vector<Term> reciever_buffer = commutator(left_term, right_term);
				mrock::utility::append_vector(reciever, std::move(reciever_buffer));
			}
		}
		return reciever;
	}

	std::ostream& operator<<(std::ostream& os, const Term& term)
	{
		if (term.multiplicity > 0) {
			os << "+";
		}
		os << term.multiplicity << " ";
		os << term.sums;
		os << term.coefficients << " ";
		for (const auto& delta : term.delta_momenta) {
			os << delta;
		}
		for (const auto& delta : term.delta_indizes) {
			os << delta;
		}
		if (term.is_identity()) {
			os << " \\hat{1} ";
			return os;
		}
		for (const auto& op : term.operators) {
			os << op << " ";
		}
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const std::vector<Term>& terms)
	{
		for (std::vector<Term>::const_iterator it = terms.begin(); it != terms.end(); ++it)
		{
			os << "\t&" << *it;
			if (it != terms.end() - 1) {
				os << " \\\\";
			}
			os << "\n";
		}
		return os;
	}

	void clean_up(std::vector<Term>& terms)
	{
		for (auto it = terms.begin(); it != terms.end();) {
			if (!(it->resolve_deltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discard_zero_momenta();
			it->rename_sums();
			it->sort();
			++it;
		}

		// Setup so that we always have a structure like delta_(l,k+something)
		for (auto& term : terms) {
			for (auto& delta : term.delta_momenta) {
				assert(delta.first.momentum_list.size() == 1U);
				int l_is = delta.first.is_used_at('l');
				if(l_is == 0) continue;

				l_is = delta.second.is_used_at('l');
				if(l_is == -1) {
					std::cout << term << std::endl;
					throw;
				}
				const Momentum l_mom('l', delta.second.momentum_list[l_is].factor);
				const Momentum remainder = delta.second - l_mom;
				delta -= remainder;
				std::swap(delta.first, delta.second);
				if(delta.first.add_Q) {
					delta.second.add_Q = !delta.second.add_Q;
					delta.first.add_Q = false;
				}
			}
		}

		clear_duplicates(terms);
		
		// Sort terms
		for (size_t i = 0; i < terms.size(); i++)
		{
			for (size_t j = i + 1; j < terms.size(); j++)
			{
				if (terms[i].sums.momenta.empty() && terms[j].sums.momenta.size() > 0) {
					std::swap(terms[i], terms[j]);
				}
				if (terms[i].sums.momenta.size() > 0 && terms[j].sums.momenta.size() > 0) {
					if (terms[i].sums.momenta.size() < terms[j].sums.momenta.size()) {
						std::swap(terms[i], terms[j]);
					}
					else if (terms[i].sums.momenta.size() == terms[j].sums.momenta.size()) {
						if (terms[i].coefficients.size() > 0) {
							if (terms[j].coefficients[0].name < terms[i].coefficients[0].name) {
								std::swap(terms[i], terms[j]);
							}
						}
					}
				}
				else if (terms[i].sums.momenta.empty() && terms[j].sums.momenta.empty()) {
					if (terms[i].coefficients.size() > 0) {
						if (terms[j].coefficients[0].name < terms[i].coefficients[0].name) {
							std::swap(terms[i], terms[j]);
						}
					}
				}
			}
		}

		for (const auto& term : terms) {
			assert(term.is_normal_ordered());
		}
	}

	void clear_duplicates(std::vector<Term>& terms) {
		// remove duplicates
		for (int i = 0; i < terms.size(); i++)
		{
			for (int j = i + 1; j < terms.size(); j++)
			{
				if (terms[i] == terms[j]) {
					terms[i].multiplicity += terms[j].multiplicity;
					terms.erase(terms.begin() + j);
					--i;
					break;
				}
			}
		}
		// removes any terms that have a 0 prefactor
		for (auto it = terms.begin(); it != terms.end();)
		{
			if (it->multiplicity == 0) {
				it = terms.erase(it);
			}
			else {
				++it;
			}
		}
	}

	std::string to_string_without_prefactor(const std::vector<Term>& terms) {
		std::string ret = "";
		for (size_t i = 0U; i < terms.size(); i++)
		{
			if (terms[i].multiplicity < 0) {
				ret += "-";
			}
			else if (i > 0) {
				ret += "+";
			}
			ret += terms[i].to_string_without_prefactor();
		}
		return ret;
	}
}
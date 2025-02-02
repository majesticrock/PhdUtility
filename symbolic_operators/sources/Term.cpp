#include <mrock/symbolic_operators/Term.hpp>
#include <mrock/symbolic_operators/KroneckerDeltaUtility.hpp>
#include <mrock/utility/RangeUtility.hpp>
#include <sstream>

namespace mrock::symbolic_operators {
	Term::Term(IntFractional _multiplicity, std::vector<Coefficient> _coefficients, const SumContainer& _sums, const std::vector<Operator>& _operators)
		: coefficients(_coefficients), sums(_sums), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const SumContainer& _sums, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), sums(_sums), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const MomentumSum& _sum_momenta, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), sums{ _sum_momenta, {} }, operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const IndexSum& _sum_indizes, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), sums{ {}, _sum_indizes }, operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, Coefficient _coefficient, const std::vector<Operator>& _operators)
		: coefficients(1, _coefficient), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, const SumContainer& _sums, const std::vector<Operator>& _operators)
		: coefficients(), sums(_sums), operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, const MomentumSum& _sum_momenta, const std::vector<Operator>& _operators)
		: coefficients(), sums{ _sum_momenta, {} }, operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, const IndexSum& _sum_indizes, const std::vector<Operator>& _operators)
		: coefficients(), sums{ {}, _sum_indizes }, operators(_operators), multiplicity(_multiplicity) {}
	Term::Term(IntFractional _multiplicity, const std::vector<Operator>& _operators)
		: coefficients(), operators(_operators), multiplicity(_multiplicity) {}

	void Term::print() const {
		std::cout << *this << std::endl;
	}

	bool Term::setDeltas()
	{
		IF_IS_TERM_TRACKED( std::cout << "setDeltas() 1:&" << (*this) << "\\\\" << std::endl; );

		for (auto& delta : delta_indizes) {
			if (!is_mutable(delta.first) && !is_mutable(delta.second)) {
				if (delta.first != delta.second) {
					return false;
				}
			}
		}

		for (auto& delta : delta_momenta)
		{
			remove_double_occurances(delta);
			if (delta.first.momentum_list.empty()) {
				if (delta.second.momentum_list.empty()) continue;
				std::swap(delta.first, delta.second);
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}
			if (delta.first.momentum_list.front().factor < 0) {
				delta.first.flip_momentum();
				delta.second.flip_momentum();
			}
			if (delta.first.momentum_list.size() > 1U && delta.second.momentum_list.empty()) {
				delta.second.momentum_list.push_back(delta.first.momentum_list[1]);
				delta.second.flip_momentum();
				delta.first.momentum_list.erase(delta.first.momentum_list.begin() + 1);
			}
		}

		IF_IS_TERM_TRACKED( std::cout << "setDeltas() 2:&" << (*this) << "\\\\" << std::endl; );
		// Removes delta_{0,Q} and delta_{0,0}
		for (auto it = delta_momenta.begin(); it != delta_momenta.end(); )
		{
			if (it->first.momentum_list.empty() && it->second.momentum_list.empty()) {
				// 0 = Q can never be achieved
				if (it->first.add_Q != it->second.add_Q) return false;
				it = delta_momenta.erase(it);
			}
			else {
				++it;
			}
		}

		// Set all deltas up to the same notation
		for (auto& delta : delta_momenta) {
			IF_IS_TERM_TRACKED( std::cout << "setDeltas() 3a:&" << (*this) << "\\\\" << std::endl; );
			for (auto& delta2 : delta_momenta) {
				remove_double_occurances(delta2);
			}
			IF_IS_TERM_TRACKED( std::cout << "setDeltas() 3b:&" << (*this) << "\\\\" << std::endl; );
			// Make sure that the first entry of each delta is not empty
			if (delta.first.momentum_list.empty()) {
				if (delta.second.momentum_list.empty()) continue;
				if (delta.second.momentum_list.size() == 1) {
					std::swap(delta.first, delta.second);
				}
				else {
					delta.first.momentum_list.push_back(delta.second.momentum_list.back());
					if (delta.first.momentum_list.front().factor > 0) {
						delta.second.flip_momentum();
					}
					else {
						delta.first.flip_momentum();
					}
					delta.second.momentum_list.pop_back();
				}
			}
			IF_IS_TERM_TRACKED( std::cout << "setDeltas() 3c:&" << (*this) << "\\\\" << std::endl; );

			// Make sure that the first entry of each delta is of size 1
			if (delta.second.momentum_list.size() == 1 && delta.first.momentum_list.size() > 1) {
				std::swap(delta.first, delta.second);
			}
			else if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() > 1) {
				bool foundCandidate = false;
				int index = 0;
				// Create a delta_{0, something} situation
				delta.second -= delta.first;
				delta.first.momentum_list.clear();
				delta.first.add_Q = false;

				// See, whether we can find a sum index within our delta
				for (auto m : sums.momenta)
				{
					index = delta.second.is_used_at(m);
					if (index > -1) {
						foundCandidate = true;
						if (abs(delta.second.momentum_list[index].factor) == 1) {
							break;
						}
					}
				}
				// If we could not find any, just use the first one we see
				if (!foundCandidate) index = 0;

				if (delta.second.momentum_list[index].factor > 0) {
					delta.second.flip_momentum();
				}
				delta.first.momentum_list.push_back(delta.second.momentum_list[index]);
				delta.first.flip_momentum();
				if (abs(delta.first.momentum_list[0].factor) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
				delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
			}
			
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}
			if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].factor < 0) {
				delta.first.flip_momentum();
				delta.second.flip_momentum();
			}
			IF_IS_TERM_TRACKED( std::cout << "setDeltas() 3d:&" << (*this) << "\\\\" << std::endl; );
			if (abs(delta.first.momentum_list[0].factor) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
			for (auto& op : operators) {
				op.momentum.replace_occurances(delta.first.momentum_list[0].name, delta.second);
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(delta.first.momentum_list[0].name, delta.second);
			}
			for (auto& delta2 : delta_momenta) {
				if (delta2 == delta) continue;
				delta2.first.replace_occurances(delta.first.momentum_list[0].name, delta.second);
				delta2.second.replace_occurances(delta.first.momentum_list[0].name, delta.second);
			}
			IF_IS_TERM_TRACKED( std::cout << "setDeltas() 3e:&" << (*this) << "\\\\" << std::endl; );
		}

		IF_IS_TERM_TRACKED( std::cout << "setDeltas() 4:&" << (*this) << "\\\\" << std::endl; );
		for (auto& delta : delta_indizes) {
			for (auto& op : operators) {
				for (auto it = op.indizes.begin(); it != op.indizes.end(); ++it)
				{
					if (delta.first == Index::SpinUp || delta.first == Index::SpinDown) {
						if (*it == delta.second) {
							*it = delta.first;
						}
					}
					else {
						if (*it == delta.first) {
							*it = delta.second;
						}
					}
				}
			}
		}

		// Remove delta^2
		remove_delta_squared(this->delta_indizes);
		remove_delta_squared(this->delta_momenta);

		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);
		remove_delta_is_one(this->delta_momenta);
		return true;
	}

	bool Term::computeSums() {
		auto changeAllIndizes = [&](const Index replaceWhat, const Index replaceWith) {
			for (auto& op : operators) {
				for (auto it = op.indizes.begin(); it != op.indizes.end(); ++it)
				{
					if (*it == replaceWhat) {
						*it = replaceWith;
					}
				}
			}
			for (auto& coeff : coefficients) {
				for (auto it = coeff.indizes.begin(); it != coeff.indizes.end(); ++it)
				{
					if (*it == replaceWhat) {
						*it = replaceWith;
					}
				}
			}
			};

		for (int i = 0; i < sums.spins.size(); i++)
		{
			for (int j = 0; j < delta_indizes.size(); j++)
			{
				if (delta_indizes[j].first == sums.spins[i]) {
					changeAllIndizes(sums.spins[i], delta_indizes[j].second);
					sums.spins.erase(sums.spins.begin() + i);
					delta_indizes.erase(delta_indizes.begin() + j);
					--i;
					break;
				}
				else if (delta_indizes[j].second == sums.spins[i]) {
					changeAllIndizes(sums.spins[i], delta_indizes[j].first);
					sums.spins.erase(sums.spins.begin() + i);
					delta_indizes.erase(delta_indizes.begin() + j);
					--i;
					break;
				}
			}
		}
		// Copy on purpose as to avoid racing conditions
		auto changeAllMomenta = [&](const MomentumSymbol::name_type replaceWhat, const Momentum replaceWith) {
			for (auto& op : operators) {
				op.momentum.replace_occurances(replaceWhat, replaceWith);
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(replaceWhat, replaceWith);
			}
			for (auto& delta : delta_momenta) {
				delta.first.replace_occurances(replaceWhat, replaceWith);
				delta.second.replace_occurances(replaceWhat, replaceWith);
			}
			};

		for (int i = 0; i < sums.momenta.size(); i++)
		{
			for (auto delta_it = delta_momenta.begin(); delta_it != delta_momenta.end(); ++delta_it)
			{
				int idx = delta_it->second.is_used_at(sums.momenta[i]);
				if(idx > -1) {
					std::swap(delta_it->first, delta_it->second);
				}
				else {
					idx = delta_it->first.is_used_at(sums.momenta[i]);
				}
				if(idx > -1) {
					assert(abs(delta_it->first.momentum_list[idx].factor) == 1);
					if(delta_it->first.momentum_list[idx].factor < 0) {
						delta_it->first.flip_momentum();
						delta_it->second.flip_momentum();
					}
					Momentum remainder = delta_it->first;
					remainder.momentum_list.erase(remainder.momentum_list.begin() + idx);
					(*delta_it) -= remainder;
					changeAllMomenta(sums.momenta[i], delta_it->second);

					sums.momenta.erase(sums.momenta.begin() + i);
					delta_momenta.erase(delta_it);
					--i;
					if (!(setDeltas())) return false;
					break;
				}
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
					flipSign();
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

	void Term::renameSums()
	{
		constexpr int N_BUFFER = 11;
		constexpr MomentumSymbol::name_type name_list[N_BUFFER]   = { 'q', 'p', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z' };
		constexpr MomentumSymbol::name_type buffer_list[N_BUFFER] = { ':', ';', '|', '?', '!', '.', '-', '_', '+', '/', '=' };
		for (size_t i = 0U; i < sums.momenta.size(); ++i)
		{
			if (i >= N_BUFFER) {
				std::cerr << "More than " << N_BUFFER << "momenta, time to implement this..." << std::endl;
				break;
			}
			if (sums.momenta[i] == name_list[i]) continue;

			for (auto& op : operators) {
				op.momentum.replace_occurances(sums.momenta[i], Momentum(buffer_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(sums.momenta[i], Momentum(buffer_list[i]));
			}
			sums.momenta[i] = name_list[i];
		}

		for (size_t i = 0U; i < sums.momenta.size(); ++i)
		{
			for (auto& op : operators) {
				op.momentum.replace_occurances(buffer_list[i], Momentum(name_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(buffer_list[i], Momentum(name_list[i]));
			}
		}

		if (sums.spins.size() == 1U && sums.spins.front() == Index::SigmaPrime) {
			sums.spins.front() = Index::Sigma;
			for (auto& op : operators) {
				for (auto& index : op.indizes) {
					if (index == Index::SigmaPrime) index = Index::Sigma;
				}
			}
			for (auto& coeff : coefficients) {
				for (auto& index : coeff.indizes) {
					if (index == Index::SigmaPrime) index = Index::Sigma;
				}
			}
		}
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

	void Term::invert_momentum(const MomentumSymbol::name_type what) {
		for (auto& coeff : coefficients) {
			coeff.invert_momentum(what);
		}
		for (auto& op : operators) {
			op.momentum.flip_single(what);
		}
	}

	void Term::invert_momentum_sum(const MomentumSymbol::name_type what) {
		if (std::find(sums.momenta.begin(), sums.momenta.end(), what) == sums.momenta.end()) {
			throw std::invalid_argument("You are trying to perform a sum transformation on a momentum that is not being summed over!");
		}
		invert_momentum(what);
	}

	void Term::remove_momentum_contribution(const MomentumSymbol::name_type value) {
		for (auto& coeff : coefficients) {
			coeff.remove_momentum_contribution(value);
		}
		for (auto& op : operators) {
			op.remove_momentum_contribution(value);
		}
		for (auto& delta : delta_momenta) {
			delta.first.remove_contribution(value);
			delta.second.remove_contribution(value);
		}
		std::erase_if(sums.momenta.summations, [&](const MomentumSymbol::name_type sum_idx) { return sum_idx == value; });
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
							terms[t].flipSign();
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
								make_delta(new_term.operators[i - 1].first_index(), new_term.operators[i].first_index()));
						}
						for (int c = 1; c < new_term.operators[i - 1].indizes.size(); c++)
						{
							// if the indizes are not the same we emplace a delta
							// otherwise no action is required
							if (new_term.operators[i - 1].indizes[c] != new_term.operators[i].indizes[c]) {
								other_deltas = true;
								new_term.delta_indizes.push_back(
									make_delta(new_term.operators[i - 1].indizes[c], new_term.operators[i].indizes[c]));
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
						if (other_deltas) terms.push_back(new_term);
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
	void commutator(std::vector<Term>& reciever, const Term& left, const Term& right)
	{
		reciever.resize(2);
		reciever[0] = left;
		reciever[0].multiplicity *= right.multiplicity;
		mrock::utility::append_vector(reciever[0].operators, right.operators);
		reciever[1] = right;
		reciever[1].multiplicity *= left.multiplicity;
		mrock::utility::append_vector(reciever[1].operators, left.operators);
		reciever[1].flipSign();

		fill_reciever(coefficients);
		fill_reciever(sums.momenta);
		fill_reciever(sums.spins);
		fill_reciever(delta_momenta);
		fill_reciever(delta_indizes);

		normal_order(reciever);
	}
	void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right)
	{
		reciever.reserve(2 * left.size() * right.size());
		
		for (const auto& left_term : left)
		{
			for (const auto& right_term : right)
			{
				std::vector<Term> reciever_buffer(2);
				commutator(reciever_buffer, left_term, right_term);
				mrock::utility::append_vector(reciever, std::move(reciever_buffer));
			}
		}
	}

	std::ostream& operator<<(std::ostream& os, const Term& term)
	{
		if (term.multiplicity > 0) {
			os << "+";
		}
		os << term.multiplicity << " \\cdot ";
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
#ifdef _TRACK_TERM
		CLEAR_TRACKED(terms);
		terms.back().is_tracked = true;
#endif
		for (std::vector<Term>::iterator it = terms.begin(); it != terms.end();) {
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discard_zero_momenta();
			if (!(it->computeSums())) {
				it = terms.erase(it);
				continue;
			}

			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discard_zero_momenta();
			it->renameSums();
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
}
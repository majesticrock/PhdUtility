#include "Term.hpp"
#include "KroneckerDeltaUtility.hpp"
#include <Utility/RangeUtility.hpp>
#include <sstream>

namespace SymbolicOperators {
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
		for (auto& delta : delta_momenta)
		{
			for (auto it = delta.first.momentum_list.begin(); it != delta.first.momentum_list.end(); )
			{
				int index = delta.second.isUsed(it->second);
				if (index < 0) {
					++it;
					continue;
				}

				int remainder = delta.second.momentum_list[index].first - it->first;
				if (remainder == 0) {
					delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
					it = delta.first.momentum_list.erase(it);
					continue;
				}

				delta.second.momentum_list[index].first = remainder;
				it = delta.first.momentum_list.erase(it);
			}
			if (delta.first.momentum_list.size() == 0) {
				if (delta.second.momentum_list.size() == 0) continue;
				std::swap(delta.first, delta.second);
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}
			if (delta.first.momentum_list.front().first < 0) {
				delta.first.flipMomentum();
				delta.second.flipMomentum();
			}
			if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() == 0) {
				delta.second.momentum_list.push_back(delta.first.momentum_list[1]);
				delta.second.flipMomentum();
				delta.first.momentum_list.erase(delta.first.momentum_list.begin() + 1);
			}
		}

		for (auto it = delta_momenta.begin(); it != delta_momenta.end(); )
		{
			if (it->first.momentum_list.size() == 0 && it->second.momentum_list.size() == 0) {
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
			for (auto& delta2 : delta_momenta) {
				for (auto it = delta2.first.momentum_list.begin(); it != delta2.first.momentum_list.end();) {
					int pos = delta2.second.isUsed(it->second);
					if (pos < 0) { ++it; continue; }
					it->first -= delta2.second.momentum_list[pos].first;
					if (it->first == 0) {
						it = delta2.first.momentum_list.erase(it);
						delta2.second.momentum_list.erase(delta2.second.momentum_list.begin() + pos);
						continue;
					}
					++it;
				}
			}

			if (delta.first.momentum_list.size() == 0) {
				if (delta.second.momentum_list.size() == 0) continue;
				if (delta.second.momentum_list.size() == 1) {
					std::swap(delta.first, delta.second);
				}
				else {
					delta.first.momentum_list.push_back(delta.second.momentum_list.back());
					if (delta.first.momentum_list.front().first > 0) {
						delta.second.flipMomentum();
					}
					else {
						delta.first.flipMomentum();
					}
					delta.second.momentum_list.pop_back();
				}
			}
			if (delta.second.momentum_list.size() == 1 && delta.first.momentum_list.size() > 1) {
				std::swap(delta.first, delta.second);
			}
			if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.size() > 1) {
				bool foundCandidate = false;
				int index = 0;
				delta.second -= delta.first;
				delta.first.momentum_list.clear();

				for (auto m : sums.momenta)
				{
					index = delta.second.isUsed(m);
					if (index >= 0) {
						foundCandidate = true;
						if (abs(delta.second.momentum_list[index].first) == 1) {
							break;
						}
					}
				}
				if (!foundCandidate) index = 0;

				if (delta.second.momentum_list[index].first > 0) {
					delta.second.flipMomentum();
				}
				delta.first.momentum_list.push_back(delta.second.momentum_list[index]);
				delta.first.flipMomentum();
				if (abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
				delta.second.momentum_list.erase(delta.second.momentum_list.begin() + index);
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}
			if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].first < 0) {
				delta.first.flipMomentum();
				delta.second.flipMomentum();
			}

			if (abs(delta.first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta.first << std::endl;
			for (auto& op : operators) {
				op.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
			}
			for (auto& delta2 : delta_momenta) {
				if (delta2 == delta) continue;
				delta2.first.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
				delta2.second.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
			}
		}
		for (auto& delta : delta_indizes) {
			for (auto& op : operators) {
				for (auto it = op.indizes.begin(); it != op.indizes.end(); ++it)
				{
					if (delta.first == SpinUp || delta.first == SpinDown) {
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

		auto changeAllMomenta = [&](const char replaceWhat, const Momentum& replaceWith) {
			for (auto& op : operators) {
				op.momentum.replaceOccurances(replaceWhat, replaceWith);
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replaceOccurances(replaceWhat, replaceWith);
			}
			for (auto it = delta_momenta.begin(); it != delta_momenta.end();) {
				it->first.replaceOccurances(replaceWhat, replaceWith);
				it->second.replaceOccurances(replaceWhat, replaceWith);
				++it;
			}
			};

		for (int i = 0; i < sums.momenta.size(); i++)
		{
			for (int j = 0; j < delta_momenta.size(); j++)
			{
				if (delta_momenta[j].first.momentum_list[0].second == sums.momenta[i]) {
					changeAllMomenta(sums.momenta[i], delta_momenta[j].second);
					if (!(delta_momenta[j].first.momentum_list.empty())) {
						if (abs(delta_momenta[j].first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta_momenta[j].first << std::endl;
					}

					sums.momenta.erase(sums.momenta.begin() + i);
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					if (!(setDeltas())) return false;
					break;
				}
				else {
					int index = delta_momenta[j].second.isUsed(sums.momenta[i]);
					if (index < 0) continue;

					Momentum buffer(delta_momenta[j].second.momentum_list[index].second, delta_momenta[j].second.momentum_list[index].first);
					if (abs(buffer.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << buffer << std::endl;
					delta_momenta[j].second.momentum_list.erase(delta_momenta[j].second.momentum_list.begin() + index);
					delta_momenta[j].second -= delta_momenta[j].first;

					if (buffer.momentum_list[0].first > 0) {
						delta_momenta[j].second.flipMomentum();
					}
					changeAllMomenta(sums.momenta[i], delta_momenta[j].second);

					sums.momenta.erase(sums.momenta.begin() + i);
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					if (!(setDeltas())) return false;
					break;
				}
			}
		}
		return true;
	}

	void Term::discardZeroMomenta() {
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

				if (coeff.translationalInvariance && !momentum.momentum_list.empty()) {
					if (momentum.momentum_list[0].first < 0) {
						momentum.flipMomentum();
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
				if (operators[i].isDaggered != operators[i - 1].isDaggered) continue;
				if (operators[i].isDaggered) {
					// c^+ c^+
					if (operators[i].indizes[0] == SpinUp && operators[i - 1].indizes[0] != SpinUp) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
						new_n = i;
					}
					else if (operators[i - 1].indizes[0] == SpinDown && operators[i].indizes[0] != SpinDown) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
						new_n = i;
					}
				}
				else {
					// c c
					if (operators[i].indizes[0] == SpinDown && operators[i - 1].indizes[0] != SpinDown) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
						new_n = i;
					}
					else if (operators[i - 1].indizes[0] == SpinUp && operators[i].indizes[0] != SpinUp) {
						std::swap(operators[i], operators[i - 1]);
						flipSign();
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
				if (operators[i].isDaggered != operators[i - 1].isDaggered) continue;
				if (operators[i].indizes[0] != operators[i - 1].indizes[0]) continue;

				if (operators[i - 1].momentum.momentum_list[0].second > operators[i].momentum.momentum_list[0].second) {
					std::swap(operators[i], operators[i - 1]);
					flipSign();
					new_n = i;
				}
			}
			n = new_n;
		}

		// check whether we can swap the sign of each momentum in the coefficients
		// 26.04.2024, I have no idea what I did here, nor do I know why I did what I did
		for (const auto& coeff : coefficients) {
			if (!(coeff.translationalInvariance)) return;
			if (std::any_of(coeff.momenta.begin(), coeff.momenta.end(), [](Momentum const& momentum) {
				return momentum.momentum_list.size() > 1U;
				})) return;
		}

		for (const auto& sum_mom : sums.momenta) {
			bool first_occurance = true;
			for (auto& op : operators) {
				int i = op.momentum.isUsed(sum_mom);
				if (i > -1) {
					if (first_occurance) {
						if (op.momentum.momentum_list[i].first < 0) {
							first_occurance = false;
						}
						else {
							break;
						}
					}
					op.momentum.momentum_list[i].first *= -1;
				}
			}
		}
	}

	void Term::renameSums()
	{
		constexpr char name_list[3] = { 'q', 'p', 'r' };
		constexpr char buffer_list[3] = { ':', ';', '|' };
		for (size_t i = 0U; i < sums.momenta.size(); ++i)
		{
			if (i >= 3) {
				std::cerr << "More than 3 momenta, time to implement this..." << std::endl;
				break;
			}
			if (sums.momenta[i] == name_list[i]) continue;

			for (auto& op : operators) {
				op.momentum.replaceOccurances(sums.momenta[i], Momentum(buffer_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replaceOccurances(sums.momenta[i], Momentum(buffer_list[i]));
			}
			sums.momenta[i] = name_list[i];
		}

		for (size_t i = 0U; i < sums.momenta.size(); ++i)
		{
			for (auto& op : operators) {
				op.momentum.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
			}
		}

		if (sums.spins.size() == 1U && sums.spins.front() == SigmaPrime) {
			sums.spins.front() = Sigma;
			for (auto& op : operators) {
				for (auto& index : op.indizes) {
					if (index == SigmaPrime) index = Sigma;
				}
			}
			for (auto& coeff : coefficients) {
				for (auto& index : coeff.indizes) {
					if (index == SigmaPrime) index = Sigma;
				}
			}
		}
	}

	std::string Term::toStringWithoutPrefactor() const
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

		if (this->isIdentity()) {
			os << " \\mathbb{1} ";
			return os.str();
		}
		for (const auto& op : this->operators) {
			os << op << " ";
		}
		return os.str();
	}

	void normalOrder(std::vector<Term>& terms) {
		for (int t = 0; t < terms.size();) {
		normalOder_outerLoop:
			if (t >= terms.size()) break;
			size_t n = terms[t].operators.size();
			size_t new_n;
			while (n > 1U) {
				new_n = 0U;
				for (size_t i = 1U; i < terms[t].operators.size(); ++i)
				{
					if (!(terms[t].operators[i - 1].isDaggered) && (terms[t].operators[i].isDaggered)) {
						bool other_deltas = false;
						new_n = i;
						// Swap cc^+
						terms[t].flipSign();
						std::swap(terms[t].operators[i - 1], terms[t].operators[i]);

						// Add a new term where cc^+ is replaced by the appropriate delta
						Term new_term(terms[t]);
						new_term.flipSign();
						if (new_term.operators[i - 1].indizes.size() != new_term.operators[i].indizes.size()) {
							throw std::invalid_argument("Operators do not have the same index count.");
						}

						if ((new_term.operators[i - 1].indizes[0] == SpinUp || new_term.operators[i - 1].indizes[0] == SpinDown)
							&& (new_term.operators[i].indizes[0] == SpinUp || new_term.operators[i].indizes[0] == SpinDown)) {
							if (new_term.operators[i - 1].indizes[0] != new_term.operators[i].indizes[0]) {
								// Case: one up the other down, then free anticommutation
								continue;
							}
						}
						else {
							new_term.delta_indizes.push_back(
								make_delta(new_term.operators[i - 1].indizes[0], new_term.operators[i].indizes[0]));
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
						// two identical fermion operators = 0
						terms.erase(terms.begin() + t);
						goto normalOder_outerLoop;
					}
				}
				n = new_n;
			}
			++t;
		}
	}

#define fill_reciever(x) reciever[0].x = left.x; Utility::append_vector(reciever[0].x, right.x); reciever[1].x = left.x; Utility::append_vector(reciever[1].x, right.x);
	void commutator(std::vector<Term>& reciever, const Term& left, const Term& right)
	{
		reciever.resize(2);
		reciever[0] = left;
		reciever[0].multiplicity *= right.multiplicity;
		Utility::append_vector(reciever[0].operators, right.operators);
		reciever[1] = right;
		reciever[1].multiplicity *= left.multiplicity;
		Utility::append_vector(reciever[1].operators, left.operators);
		reciever[1].flipSign();

		fill_reciever(coefficients);
		fill_reciever(sums.momenta);
		fill_reciever(sums.spins);
		fill_reciever(delta_momenta);
		fill_reciever(delta_indizes);

		normalOrder(reciever);
	}
	void commutator(std::vector<Term>& reciever, const std::vector<Term>& left, const std::vector<Term>& right)
	{
		reciever.reserve(2 * left.size() * right.size());
		std::vector<Term> reciever_buffer(2);

		for (int i = 0; i < left.size(); i++)
		{
			for (int j = 0; j < right.size(); j++)
			{
				commutator(reciever_buffer, left[i], right[j]);
				Utility::append_vector(reciever, reciever_buffer);
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
		if (term.isIdentity()) {
			os << " \\mathbb{1} ";
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

	void cleanUp(std::vector<Term>& terms)
	{
		for (std::vector<Term>::iterator it = terms.begin(); it != terms.end();) {
			//std::cout << count++ << " of " << terms.size() << ":&\t" << *it << "\\\\" << std::endl;
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			if (!(it->computeSums())) {
				it = terms.erase(it);
				continue;
			}
			if (!(it->setDeltas())) {
				it = terms.erase(it);
				continue;
			}
			it->discardZeroMomenta();
			it->renameSums();
			it->sort();

			++it;
		}

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
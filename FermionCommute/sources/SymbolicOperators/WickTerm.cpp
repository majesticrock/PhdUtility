#include "WickTerm.hpp"
#include <Utility/StringUtility.hpp>
#include "KroneckerDeltaUtility.hpp"
#include <cctype>
#include <cassert>

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].indizes[0]
#define R_SPIN temporary_operators[i + 1].indizes[0]

namespace SymbolicOperators {
	// Constructors
	WickTerm::WickTerm(const Term* base)
		: multiplicity(base->multiplicity), coefficients(base->coefficients), sums(base->sums), operators(),
		delta_momenta(base->delta_momenta), delta_indizes(base->delta_indizes), temporary_operators()
	{
	}
	WickTerm::WickTerm(const Term& base)
		: multiplicity(base.multiplicity), coefficients(base.coefficients), sums(base.sums), operators(),
		delta_momenta(base.delta_momenta), delta_indizes(base.delta_indizes), temporary_operators()
	{
	}
	WickTerm::WickTerm(const WickTerm& base, const TemplateResult::SingleResult& result)
		: multiplicity(result.factor* base.multiplicity), coefficients(base.coefficients), sums(base.sums), operators(base.operators),
		delta_momenta(base.delta_momenta), delta_indizes(base.delta_indizes), temporary_operators()
	{
		this->operators.push_back(result.op);
		this->delta_indizes.insert(this->delta_indizes.end(), result.index_deltas.begin(), result.index_deltas.end());
	}
	WickTerm::WickTerm(const std::string& expression) : multiplicity(1)
	{
		// Syntax
		// [factor] [index_sum] [momentum_sum] [coefficients...] [momentum_deltas...] [index_deltas...] [operators...]
		/* factor needs to be an integer
		*  Sums must be "sum:index{index1,index2,...}" or "sum:momentum{momentum_name1,momentum_name2,...}"
		*  coefficient must be "c:name{Momentum_expression1,...;index1,index2,...}"
		*  deltas must be "delta:momentum{Momentum_expression,Momentum_expression}" or "delta:index{Index,Index}"
		*  operators must be "o:type{Momentum_expression;index1,index2,...}(^+)"
		*/

		size_t pos{};
		if (std::isdigit(expression.front()) || expression.front() == '-' || expression.front() == '+') {
			pos = expression.find(' ');
			if (pos != std::string::npos) {
				this->multiplicity = std::stoi(expression.substr(0U, pos));
			}
		}
		size_t new_pos{};
		++pos;
		while (new_pos != std::string::npos) {
			new_pos = expression.find(' ', pos);
			string_parser(expression.substr(pos, new_pos - pos));
			pos = new_pos + 1;
		}
	}

	// Member functions
	void WickTerm::string_parser(std::string&& expression) {
		size_t forward_pos{ expression.find(' ') };
		if (forward_pos == std::string::npos)
			forward_pos = expression.size();

		const std::string sub = expression.substr(0U, forward_pos);
		const size_t sub_delimiter = sub.find(':');

		if (sub_delimiter == std::string::npos)
			throw std::invalid_argument("Did not find ':' in " + expression);

		if (sub.substr(0U, sub_delimiter) == "sum") {
			const std::string type = expression.substr(sub_delimiter + 1, expression.find('{', sub_delimiter) - sub_delimiter - 1);
			const std::vector<std::string> argument_list = Utility::extract_elements(expression);

			if (type == "index") {
				this->sums.spins.reserve(argument_list.size());
				for (const auto& arg : argument_list) {
					this->sums.spins.push_back(string_to_index.at(arg));
				}
			}
			else if (type == "momentum") {
				this->sums.momenta.reserve(argument_list.size());
				for (const auto& arg : argument_list) {
					assert(arg.size() == 1U);
					this->sums.momenta.push_back(arg.front());
				}
			}
			else {
				throw std::invalid_argument("Sum type not recognized " + type + " in expression " + expression);
			}
		}
		else if (sub.substr(0U, sub_delimiter) == "delta") {
			const std::string type = expression.substr(sub_delimiter + 1, expression.find('{', sub_delimiter) - sub_delimiter - 1);
			const std::vector<std::string> argument_list = Utility::extract_elements(expression);
			assert(argument_list.size() == 2U);

			if (type == "index") {
				this->delta_indizes.push_back(make_delta(string_to_index.at(argument_list[0]), string_to_index.at(argument_list[1])));
			}
			else if (type == "momentum") {
				this->delta_momenta.push_back(make_delta(Momentum(argument_list[0]), Momentum(argument_list[1])));
			}
			else {
				throw std::invalid_argument("Delta type not recognized " + type + " in expression " + expression);
			}
		}
		else if (sub.substr(0U, sub_delimiter) == "c") {
			this->coefficients.push_back(Coefficient::parse_string(sub.substr(sub_delimiter + 1)));
		}
		else if (sub.substr(0U, sub_delimiter) == "o") {
			this->operators.push_back(WickOperator(sub.substr(sub_delimiter + 1)));
		}
		else {
			throw std::invalid_argument("Did parse expression <" + expression + "> at <" + sub + "> with delimiter " + std::to_string(sub_delimiter));
		}
	}

	bool WickTerm::setDeltas()
	{
		//remove_delta_squared(this->delta_indizes);
		//remove_delta_squared(this->delta_momenta);

		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);
		remove_delta_is_one(this->delta_momenta);

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
			if (delta.first.momentum_list.empty()) {
				if (delta.second.momentum_list.empty()) continue;
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
			if (delta.first.momentum_list.size() > 1 && delta.second.momentum_list.empty()) {
				delta.second.momentum_list.push_back(delta.first.momentum_list[1]);
				delta.second.flipMomentum();
				delta.first.momentum_list.erase(delta.first.momentum_list.begin() + 1);
			}
		}

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

			if (delta.first.momentum_list.empty()) {
				if (delta.second.momentum_list.empty()) continue;
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

			if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].first < 0) {
				delta.first.flipMomentum();
				delta.second.flipMomentum();
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
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

		//remove_delta_squared(this->delta_indizes);
		//remove_delta_squared(this->delta_momenta);
		// Remove delta^2
		for (int i = 0; i < delta_momenta.size(); i++)
		{
			for (int j = i + 1; j < delta_momenta.size(); j++)
			{
				if (delta_momenta[i] == delta_momenta[j]) {
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					break;
				}

				auto delta_buffer = delta_momenta[j];
				delta_buffer.first.flipMomentum();
				delta_buffer.second.flipMomentum();
				if (delta_momenta[i] == delta_buffer) {
					delta_momenta.erase(delta_momenta.begin() + j);
					--i;
					break;
				}
			}
		}
		for (int i = 0; i < delta_indizes.size(); i++)
		{
			for (int j = i + 1; j < delta_indizes.size(); j++)
			{
				if (delta_indizes[i] == delta_indizes[j]) {
					delta_indizes.erase(delta_indizes.begin() + j);
					--i;
					break;
				}
			}
		}

		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);
		remove_delta_is_one(this->delta_momenta);

		return !(is_always_zero(this->delta_indizes) || is_always_zero(this->delta_momenta));
	}

	bool WickTerm::computeSums()
	{
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
			for (auto& delta : delta_indizes) {
				if (delta.first == replaceWhat) {
					delta.first = replaceWith;
				}
				if (delta.second == replaceWhat) {
					delta.second = replaceWith;
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

		auto changeAllMomenta = [&](const char replaceWhat, const Momentum replaceWith) {
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
					if (abs(delta_momenta[j].first.momentum_list[0].first) != 1) std::cerr << "Not yet implemented! " << delta_momenta[j].first << std::endl;
					changeAllMomenta(sums.momenta[i], delta_momenta[j].second);

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

	void WickTerm::discardZeroMomenta()
	{
		for (auto& op : operators) {
			op.momentum.remove_zeros();
		}
		for (auto& coeff : coefficients) {
			coeff.momenta.remove_zeros();
		}
	}

	void WickTerm::renameSums()
	{
		constexpr char name_list[3] = { 'q', 'p', 'r' };
		constexpr char buffer_list[3] = { ':', ';', '|' };
		for (int i = 0; i < sums.momenta.size(); i++)
		{
			if (i >= 3) {
				throw std::invalid_argument("More than 3 momenta, time to implement this...");
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

		for (int i = 0; i < sums.momenta.size(); i++)
		{
			for (auto& op : operators) {
				op.momentum.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replaceOccurances(buffer_list[i], Momentum(name_list[i]));
			}
		}

		for (const auto& sum : sums.momenta)
		{
			for (auto& op : operators) {
				int index = op.momentum.isUsed(sum);
				if (index < 0) continue;
				if (op.momentum.momentum_list.size() == 1) break;

				Momentum buffer = op.momentum;
				if (buffer.momentum_list[index].first > 0) buffer.flipMomentum();
				buffer.momentum_list[index].first *= -1;
				buffer.momentum_list[index].second = buffer_list[0];

				for (auto& op2 : operators) {
					op2.momentum.replaceOccurances(sum, buffer);
					op2.momentum.replaceOccurances(buffer_list[0], Momentum(sum));
				}
				for (auto& coeff : coefficients) {
					coeff.momenta.replaceOccurances(sum, buffer);
					coeff.momenta.replaceOccurances(buffer_list[0], Momentum(sum));
				}
			}
		}
		discardZeroMomenta();

		if(sums.spins.size() == 1U && sums.spins.front() == SigmaPrime) {
			for(auto& coeff : coefficients){
				for(auto& index : coeff.indizes){
					if(index == SigmaPrime) index = Sigma;
				}
			}
			for(auto& op : operators){
				for(auto& index : op.indizes){
					if(index == SigmaPrime) index = Sigma;
				}
			}
			for(auto& delta_index : delta_indizes){
				if(delta_index.first == SigmaPrime) delta_index.first = Sigma;
				if(delta_index.second == SigmaPrime) delta_index.second = Sigma;
			}
			sums.spins.front() = Sigma;
		}
	}

	void WickTerm::sort()
	{
		for (auto& delta : delta_momenta) {
			if (delta.first.momentum_list.size() == 1 && delta.second.momentum_list.size() == 1) {
				// This comparison is well defined because we save the momentum as char i.e. byte
				// which is easily comparable
				if (delta.first.momentum_list[0].second < delta.second.momentum_list[0].second) {
					std::swap(delta.first, delta.second);
					if (delta.first.momentum_list[0].first < 0) {
						delta.first.flipMomentum();
						delta.second.flipMomentum();
					}
					if (delta.first.add_Q) {
						delta.first.add_Q = false;
						delta.second.add_Q = !(delta.second.add_Q);
					}
				}
				for (auto& op : operators) {
					op.momentum.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
				}
				for (auto& coeff : coefficients) {
					coeff.momenta.replaceOccurances(delta.first.momentum_list[0].second, delta.second);
				}
			}
		}

		for (auto& op : operators) {
			if (op.type == CDW_Type && op.momentum.add_Q) {
				op.momentum.add_Q = false;
				op.isDaggered = !(op.isDaggered);
			}
		}

		for (int i = 0U; i < operators.size(); ++i)
		{
			for (int j = i + 1U; j < operators.size(); ++j)
			{
				if (operators[i].type > operators[j].type) {
					std::swap(operators[i], operators[j]);
				}
				else if (operators[i].type == operators[j].type) {
					if (operators[i].momentum.momentum_list[0].second > operators[j].momentum.momentum_list[0].second) {
						std::swap(operators[i], operators[j]);
					}
					else if (operators[i].momentum.momentum_list[0].second == operators[j].momentum.momentum_list[0].second) {
						if (operators[i].momentum.add_Q && !(operators[j].momentum.add_Q)) {
							std::swap(operators[i], operators[j]);
						}
					}
				}
			}
		}

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
					this->multiplicity *= -1;
				}
			}
		}

		for (auto& coeff : coefficients) {
			for (auto& momentum : coeff.momenta) {
				for (const auto& sum : sums.momenta) {
					int idx = momentum.isUsed(sum);
					if (idx < 0) continue;

					if (momentum.momentum_list[idx].first < 0) {
						for (auto& op : operators) {
							op.momentum.flip_single(sum);
						}
						for (auto& coeff2 : coefficients) {
							coeff2.momenta.flip_single(sum);
						}
					}
				}
			}
		}
	}

	void WickTerm::includeTemplateResult(const TemplateResult::SingleResult& result) {
		this->delta_indizes.insert(this->delta_indizes.begin(), result.index_deltas.begin(), result.index_deltas.end());
		this->operators.push_back(result.op);
		this->multiplicity *= result.factor;
	}
	
	// Operator overloads
	std::ostream& operator<<(std::ostream& os, const WickTerm& term)
	{
		if (term.multiplicity > 0) {
			os << "+";
		}
		os << term.multiplicity << " \\cdot ";
		os << term.sums;
		os << term.coefficients << " ";
		for (const auto& delta : term.delta_momenta) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}
		for (const auto& delta : term.delta_indizes) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
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
	std::ostream& operator<<(std::ostream& os, const WickTermCollector& terms)
	{
		for (WickTermCollector::const_iterator it = terms.begin(); it != terms.end(); ++it)
		{
			os << "\t&" << *it;
			if (it != terms.end() - 1) {
				os << " \\\\";
			}
			os << "\n";
		}
		return os;
	}
	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTerm& rhs)
	{
		for (auto it = lhs.begin(); it != lhs.end(); ++it) {
			if (*it == rhs) {
				it->multiplicity += rhs.multiplicity;
				if (it->multiplicity == 0)
					lhs.erase(it);
				return lhs;
			}
		}
		lhs.push_back(rhs);
		return lhs;
	}
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTerm& rhs)
	{
		for (auto it = lhs.begin(); it != lhs.end(); ++it) {
			if (*it == rhs) {
				it->multiplicity -= rhs.multiplicity;
				if (it->multiplicity == 0)
					lhs.erase(it);
				return lhs;
			}
		}
		lhs.push_back(rhs);
		return lhs;
	}
	WickTermCollector& operator+=(WickTermCollector& lhs, const WickTermCollector& rhs)
	{
		for (const auto& term : rhs) {
			lhs += term;
		}
		return lhs;
	}
	WickTermCollector& operator-=(WickTermCollector& lhs, const WickTermCollector& rhs)
	{
		for (const auto& term : rhs) {
			lhs -= term;
		}
		return lhs;
	}
}
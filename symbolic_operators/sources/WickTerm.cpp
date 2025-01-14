#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/symbolic_operators/KroneckerDeltaUtility.hpp>
#include <mrock/utility/StringUtility.hpp>
#include <cctype>
#include <cassert>

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].first_index()
#define R_SPIN temporary_operators[i + 1].first_index()

namespace mrock::symbolic_operators {
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
			const std::vector<std::string> argument_list = mrock::utility::extract_elements(expression);

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
			const std::vector<std::string> argument_list = mrock::utility::extract_elements(expression);
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
			throw std::invalid_argument("Did not parse expression <" + expression + "> at <" + sub + "> with delimiter " + std::to_string(sub_delimiter));
		}
	}

	bool WickTerm::setDeltas()
	{
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);
		remove_delta_is_one(this->delta_momenta);

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
			for (auto& delta2 : delta_momenta) {
				remove_double_occurances(delta2);
			}
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
					index = delta.second.isUsed(m);
					if (index >= 0) {
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

			if (delta.first.momentum_list.size() == 1 && delta.first.momentum_list[0].factor < 0) {
				delta.first.flip_momentum();
				delta.second.flip_momentum();
			}
			if (delta.first.add_Q) {
				delta.first.add_Q = false;
				delta.second.add_Q = !(delta.second.add_Q);
			}

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
		}
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
				delta_buffer.first.flip_momentum();
				delta_buffer.second.flip_momentum();
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

		auto changeAllMomenta = [&](const MomentumSymbol::name_type replaceWhat, const Momentum replaceWith) {
			for (auto& op : operators) {
				op.momentum.replace_occurances(replaceWhat, replaceWith);
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(replaceWhat, replaceWith);
			}
			for (auto it = delta_momenta.begin(); it != delta_momenta.end();) {
				it->first.replace_occurances(replaceWhat, replaceWith);
				it->second.replace_occurances(replaceWhat, replaceWith);
				++it;
			}
			};

		for (int i = 0; i < sums.momenta.size(); i++)
		{
			for (int j = 0; j < delta_momenta.size(); j++)
			{
				if (delta_momenta[j].first.momentum_list[0].name == sums.momenta[i]) {
					if (abs(delta_momenta[j].first.momentum_list[0].factor) != 1) std::cerr << "Not yet implemented! " << delta_momenta[j].first << std::endl;
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

					Momentum buffer(delta_momenta[j].second.momentum_list[index].name, delta_momenta[j].second.momentum_list[index].factor);
					if (abs(buffer.momentum_list[0].factor) != 1) std::cerr << "Not yet implemented! " << buffer << std::endl;
					delta_momenta[j].second.momentum_list.erase(delta_momenta[j].second.momentum_list.begin() + index);
					delta_momenta[j].second -= delta_momenta[j].first;

					if (buffer.momentum_list[0].factor > 0) {
						delta_momenta[j].second.flip_momentum();
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

	void WickTerm::discard_zero_momenta()
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
		constexpr MomentumSymbol::name_type name_list[3] = { 'q', 'p', 'r' };
		constexpr MomentumSymbol::name_type buffer_list[3] = { ':', ';', '|' };
		for (int i = 0; i < sums.momenta.size(); i++)
		{
			if (i >= 3) {
				throw std::invalid_argument("More than 3 momenta, time to implement this...");
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

		for (int i = 0; i < sums.momenta.size(); i++)
		{
			for (auto& op : operators) {
				op.momentum.replace_occurances(buffer_list[i], Momentum(name_list[i]));
			}
			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(buffer_list[i], Momentum(name_list[i]));
			}
		}

		for (const auto& sum : sums.momenta)
		{
			for (auto& op : operators) {
				int index = op.momentum.isUsed(sum);
				if (index < 0) continue;
				if (op.momentum.momentum_list.size() == 1) break;

				Momentum buffer = op.momentum;
				if (buffer.momentum_list[index].factor > 0) buffer.flip_momentum();
				buffer.momentum_list[index].factor *= -1;
				buffer.momentum_list[index].name = buffer_list[0];

				for (auto& op2 : operators) {
					op2.momentum.replace_occurances(sum, buffer);
					op2.momentum.replace_occurances(buffer_list[0], Momentum(sum));
				}
				for (auto& coeff : coefficients) {
					coeff.momenta.replace_occurances(sum, buffer);
					coeff.momenta.replace_occurances(buffer_list[0], Momentum(sum));
				}
			}
		}
		discard_zero_momenta();

		if(sums.spins.size() == 1U && sums.spins.front() == Index::SigmaPrime) {
			for (auto& coeff : coefficients) {
				for (auto& index : coeff.indizes) {
					if(index == Index::SigmaPrime) index = Index::Sigma;
				}
			}
			for (auto& op : operators) {
				for (auto& index : op.indizes) {
					if(index == Index::SigmaPrime) index = Index::Sigma;
				}
			}
			for (auto& delta_index : delta_indizes) {
				if(delta_index.first == Index::SigmaPrime) delta_index.first = Index::Sigma;
				if(delta_index.second == Index::SigmaPrime) delta_index.second = Index::Sigma;
			}
			sums.spins.front() = Index::Sigma;
		}
	}

	void WickTerm::sort()
	{
		for (auto& delta : delta_momenta) {
			if (delta.first.momentum_list.size() == 1 && delta.second.momentum_list.size() == 1) {
				// This comparison is well defined because we save the momentum as char i.e. byte
				// which is easily comparable
				if (delta.first.momentum_list[0].name < delta.second.momentum_list[0].name) {
					std::swap(delta.first, delta.second);
					if (delta.first.momentum_list[0].factor < 0) {
						delta.first.flip_momentum();
						delta.second.flip_momentum();
					}
					if (delta.first.add_Q) {
						delta.first.add_Q = false;
						delta.second.add_Q = !(delta.second.add_Q);
					}
				}
				for (auto& op : operators) {
					op.momentum.replace_occurances(delta.first.momentum_list[0].name, delta.second);
				}
				for (auto& coeff : coefficients) {
					coeff.momenta.replace_occurances(delta.first.momentum_list[0].name, delta.second);
				}
			}
		}

		for (auto& op : operators) {
			if (op.type == CDW_Type && op.momentum.add_Q) {
				op.momentum.add_Q = false;
				op.is_daggered = !(op.is_daggered);
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
					if (momentum_order(operators[i].momentum, operators[j].momentum)) {
						std::swap(operators[i], operators[j]);
					}
				}
			}
		}

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
					this->multiplicity *= -1;
				}
			}
		}

		for (auto& coeff : coefficients) {
			if (coeff.momenta.size() == 3U) {
				Momentum* first_momentum = &coeff.momenta.front();
				if (first_momentum != nullptr && first_momentum->empty()) {
					if (coeff.momenta.size() > 1U) first_momentum = &coeff.momenta[1];
				}
				if ((first_momentum != nullptr) && (!first_momentum->empty())) {
					if (!first_momentum->first_momentum_is('k')) {
						coeff.use_symmetric_interaction_exchange();
					}
					if (sums.momenta.empty()) {
						if (coeff.momenta.back().first_momentum_is_negative()) {
							coeff.use_symmetric_interaction_inversion();
						}
					}
					else if ((!sums.momenta.is_summed_over(first_momentum->front().name)) && first_momentum->first_momentum_is_negative()) {
						coeff.use_symmetric_interaction_inversion();
					}
				}
			}
			
			for (auto& momentum : coeff.momenta) {
				if (momentum.empty()) continue;
				if ((!sums.momenta.is_summed_over(momentum.front().name)) && momentum.front().factor < 0) {
					if (coeff.inversion_symmetry) momentum.flip_momentum();
				}

				for (const auto& sum : sums.momenta) {
					int idx = momentum.isUsed(sum);
					if (idx < 0) continue;

					if (momentum.momentum_list[idx].factor < 0) {
						invert_momentum_sum(sum);
					}
				}
			}
		}
	}

	void WickTerm::include_template_result(const TemplateResult::SingleResult& result) {
		this->delta_indizes.insert(this->delta_indizes.begin(), result.index_deltas.begin(), result.index_deltas.end());
		this->operators.push_back(result.op);
		this->multiplicity *= result.factor;
	}
	
	void WickTerm::invert_momentum(const MomentumSymbol::name_type what) {
		for (auto& coeff : coefficients) {
			coeff.invert_momentum(what);
		}
		for (auto& op : operators) {
			op.momentum.flip_single(what);
		}
	}

	void WickTerm::invert_momentum_sum(const MomentumSymbol::name_type what) {
		if (std::find(sums.momenta.begin(), sums.momenta.end(), what) == sums.momenta.end()) {
			throw std::invalid_argument("You are trying to perform a sum transformation on a momentum that is not being summed over!");
		}
		invert_momentum(what);
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
		if (term.is_identity()) {
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
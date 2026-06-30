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

	bool WickTerm::resolve_momentum_deltas() 
	{
		for (auto delta_it = delta_momenta.begin(); delta_it != delta_momenta.end(); ) {
			delta_it->first -= delta_it->second;
			delta_it->second = Momentum();

			if (delta_it->first.momentum_list.empty() && delta_it->second.momentum_list.empty()) {
				// 0 = Q can never be achieved
				if (delta_it->first.add_Q != delta_it->second.add_Q) return false;
				// delta_(0,0) = 1
				delta_it = delta_momenta.erase(delta_it);
				continue;
			}

			MomentumSymbol resolve_to{ *(delta_it->first.begin()) };
			bool found_sum{};

			for (auto sum_it = sums.momenta.begin(); sum_it != sums.momenta.end(); ++sum_it) {
				const auto found_it = std::find_if(delta_it->first.begin(), delta_it->first.end(), [&sum_it](const MomentumSymbol& symbol) {
					return symbol.name == *sum_it;
				});
				if ( found_it != delta_it->first.end()) {
					resolve_to = *found_it;
					sums.momenta.erase(sum_it);
					found_sum = true;
					break;
				}
			}

			if (resolve_to.factor > 0) {
				delta_it->first.flip_momentum();
			}
			else {
				resolve_to.factor *= -1;
			}
			delta_it->second = Momentum(resolve_to);
			delta_it->first += delta_it->second;
			

			for (MomentumSymbol& symbol : delta_it->first) {
				assert(symbol.factor % delta_it->second.front().factor == 0);
				symbol.factor /= delta_it->second.front().factor;
			}

			for (auto& coeff : coefficients) {
				coeff.momenta.replace_occurances(delta_it->second.front().name, delta_it->first);
			}
			for (auto& op : operators) {
				op.momentum.replace_occurances(delta_it->second.front().name, delta_it->first);
			}
			for (auto delta_it2 = delta_momenta.begin(); delta_it2 != delta_momenta.end(); ++delta_it2) {
				if (delta_it2 == delta_it) continue;
				delta_it2->first.replace_occurances(delta_it->second.front().name, delta_it->first);
				delta_it2->second.replace_occurances(delta_it->second.front().name, delta_it->first);
			}

			if (found_sum) {
				delta_it = delta_momenta.erase(delta_it);
			}
			else {
				++delta_it;
			}
		}

		// Make sure that delta.first has always exactly one momentum symbol (or delta_{0,0})
		for (auto& delta : delta_momenta) {
			if (delta.first.empty() && !delta.second.empty()) {
				std::swap(delta.first, delta.second);
			} 
			if (delta.first.size() > 1U) {
				const Momentum shift = delta.first - Momentum(delta.first.front());
				delta.first -= shift;
				delta.second -= shift;
			}
			if (delta.first.front().factor < 0) {
				delta.first.flip_momentum();
				delta.second.flip_momentum();
			}
		}

		// Remove delta^2
		remove_delta_squared(this->delta_indizes);
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);

		return true;
	}


	bool WickTerm::resolve_index_deltas() 
	{
		if (is_always_zero(delta_indizes)) return false;

		for (auto delta_it = delta_indizes.begin(); delta_it != delta_indizes.end(); ) {
			Index to_resolve { Index::UndefinedIndex };
			Index change_to { Index::UndefinedIndex };
			bool found_sum{};
			auto sum_it = std::find_if( sums.spins.begin(), sums.spins.end(), [&delta_it](const Index& idx) {
				return idx == delta_it->first;
			});
			if (sum_it != sums.spins.end()) {
				to_resolve = delta_it->first;
				change_to = delta_it->second;
				found_sum = true;
			}
			else {
				sum_it = std::find_if( sums.spins.begin(), sums.spins.end(), [&delta_it](const Index& idx) {
					return idx == delta_it->second;
				});
				if (sum_it != sums.spins.end()) {
					to_resolve = delta_it->second;
					change_to = delta_it->first;
					found_sum = true;
				}
			}

			if (to_resolve == Index::UndefinedIndex) {
				if (is_mutable(delta_it->first)) {
					to_resolve = delta_it->first;
					change_to = delta_it->second;
				}
				else if (is_mutable(delta_it->second)) {
					to_resolve = delta_it->second;
					change_to = delta_it->first;
				}
				else if(delta_it->first != delta_it->second) {
					// Two differing, immutable indizes can never be equal
					return false;
				}
			}

			if(delta_it->first == delta_it->second) continue;

			for (auto& op : operators) {
				op.indizes.replace_index(to_resolve, change_to);
			}
			for (auto& coeff : coefficients) {
				coeff.indizes.replace_index(to_resolve, change_to);
			}

			for (auto delta_it2 = delta_indizes.begin(); delta_it2 != delta_indizes.end(); ++delta_it2) {
				if (delta_it2 == delta_it) continue;
				if (delta_it2->first == to_resolve) {
					delta_it2->first = change_to;
				}
				if (delta_it2->second == to_resolve) {
					delta_it2->second = change_to;
				}
			}

			if (found_sum) {
				sums.spins.erase(sum_it);
				delta_it = delta_indizes.erase(delta_it);
			}
			else {
				++delta_it;
			}
		}

		// Remove delta^2
		remove_delta_squared(this->delta_indizes);
		// Erase delta_k,k etc
		remove_delta_is_one(this->delta_indizes);

		return true;
	}


	bool WickTerm::resolve_deltas()
	{
		if (!resolve_momentum_deltas()) return false;
		if (!resolve_index_deltas()) return false;
		
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

	void WickTerm::rename_sums()
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
				int index = op.momentum.is_used_at(sum);
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
					int idx = momentum.is_used_at(sum);
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

	bool WickTerm::is_pauli_forbidden() const {
		// Cannot be forbidden if there is only one WickOperators 
		if (this->operators.size() < 2U) return false; 

		for (auto it = this->operators.begin(); it != this->operators.end() - 1; ++it) {
			for (auto jt = it + 1; jt != this->operators.end(); ++jt) {
				if (*it == *jt) {
					return true;
				}
			}
		}
		return false;
	}

	// Operator overloads
	std::ostream& operator<<(std::ostream& os, const WickTerm& term)
	{
		if (term.multiplicity > 0) {
			os << "+";
		}
		os << term.multiplicity << " ";
		os << term.sums;
		os << term.coefficients << " ";
		for (const auto& delta : term.delta_momenta) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
		}
		for (const auto& delta : term.delta_indizes) {
			os << "\\delta_{" << delta.first << ", " << delta.second << "} ";
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
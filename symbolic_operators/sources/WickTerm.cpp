#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/symbolic_operators/KroneckerDeltaUtility.hpp>
#include <mrock/symbolic_operators/detail/string_helper.hpp>
#include <cctype>
#include <cassert>

#define LEFT temporary_operators[i]
#define RIGHT temporary_operators[i + 1]
#define L_SPIN temporary_operators[i].first_index()
#define R_SPIN temporary_operators[i + 1].first_index()

namespace mrock::symbolic_operators {
	// Constructors
	WickTerm::WickTerm(const Term* base)
		: AbstractTerm<WickOperator>(base->multiplicity, 
			base->coefficients, 
			base->sums, 
			base->delta_momenta, 
			base->delta_indizes, 
			std::vector<WickOperator>()),
		temporary_operators()
	{}

	WickTerm::WickTerm(const Term& base)
		: AbstractTerm<WickOperator>(base.multiplicity, 
			base.coefficients, 
			base.sums, 
			base.delta_momenta, 
			base.delta_indizes, 
			std::vector<WickOperator>()),
		temporary_operators()
	{}

	WickTerm::WickTerm(const WickTerm& base, const TemplateResult::SingleResult& result)
		: AbstractTerm<WickOperator>(base.multiplicity, 
			base.coefficients, 
			base.sums, 
			base.delta_momenta, 
			base.delta_indizes, 
			base.operators),
		temporary_operators()
	{
		this->operators.push_back(result.op);
		this->delta_indizes.insert(this->delta_indizes.end(), result.index_deltas.begin(), result.index_deltas.end());
	}

	WickTerm::WickTerm(const std::string& expression) : AbstractTerm<WickOperator>(1)
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
			const std::vector<std::string> argument_list = extract_elements(expression);

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
			const std::vector<std::string> argument_list = extract_elements(expression);
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

	bool WickTerm::resolve_deltas()
	{
		if (!resolve_momentum_deltas()) return false;
		if (!resolve_index_deltas()) return false;
		
		return true;
	}

	void WickTerm::rename_sums() {
		AbstractTerm<WickOperator>::rename_sums();

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
	
	bool WickTerm::is_pauli_forbidden() const {
		// Cannot be forbidden if there is only one WickOperator
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
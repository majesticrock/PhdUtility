#include "Coefficient.hpp"
#include <Utility/StringUtility.hpp>

namespace SymbolicOperators {
	Coefficient Coefficient::parse_string(const std::string& expression)
	{
		// Syntax:   name{Momentum_expression1,Momentum_expression1;index1,index2,...}

		Coefficient ret;
		ret.name = expression.substr(0U, expression.find('{'));
		std::vector<std::string> momentum_strings = Utility::extract_elements(expression, '{', ';');
		std::vector<std::string> index_strings = Utility::extract_elements(expression, ';', '}');

		ret.momenta.reserve(momentum_strings.size());
		for (const auto& arg : momentum_strings) {
			ret.momenta.push_back(Momentum(arg));
		}
		ret.indizes.reserve(index_strings.size());
		for (const auto& arg : index_strings) {
			ret.indizes.push_back(string_to_index.at(arg));
		}

		return ret;
	}

	std::ostream& operator<<(std::ostream& os, const Coefficient& coeff)
	{
		os << coeff.name;
		if (!coeff.indizes.empty()) {
			os << "_{ " << coeff.indizes << "}";
		}
		if (coeff.isDaggered) {
			os << "^*";
		}
		os << coeff.momenta << " ";
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs) {
		for (std::vector<Coefficient>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it)
		{
			os << (*it) << " ";
		}
		return os;
	}

	Coefficient::Coefficient()
		: name(""), momenta(), indizes(), Q_changes_sign(false), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name)
		: name(_name), momenta(), indizes(), Q_changes_sign(false), isDaggered(false) {}
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momenta(_momentum), indizes(_indizes), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, const Momentum& _momentum, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momenta(_momentum), indizes(), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) {}
	Coefficient::Coefficient(std::string _name, const MomentumList& _momenta, const IndexWrapper& _indizes, bool _Q_changes_sign, bool _isDaggered)
		: name(_name), momenta(_momenta), indizes(), Q_changes_sign(_Q_changes_sign), isDaggered(_isDaggered) { }
}
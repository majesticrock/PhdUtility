#include "WickOperator.hpp"
#include <Utility/StringUtility.hpp>
#include <cassert>

namespace SymbolicOperators {
	WickOperator::WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const IndexWrapper& _indizes)
		: type(_type), isDaggered(_isDaggered), momentum(_momentum), indizes(_indizes) {}
	WickOperator::WickOperator(const OperatorType& _type, const bool _isDaggered, const Momentum& _momentum, const Index _index)
		: type(_type), isDaggered(_isDaggered), momentum(_momentum), indizes(_index) {}
	WickOperator::WickOperator()
		: type(Undefined_Type), isDaggered(false), momentum(), indizes() {}

	WickOperator::WickOperator(const std::string& expression)
	{
		// Syntax    type{Momentum_expression;index1,index2,...}(^+)

		this->type = string_to_wick.at(expression.substr(0U, expression.find('{')));
		std::vector<std::string> momentum_strings = Utility::extract_elements(expression, '{', ';');
		std::vector<std::string> index_strings = Utility::extract_elements(expression, ';', '}');

		assert(momentum_strings.size() == 1U);
		this->momentum = Momentum(momentum_strings.front());

		this->indizes.reserve(index_strings.size());
		for (const auto& arg : index_strings) {
			this->indizes.push_back(string_to_index.at(arg));
		}

		this->isDaggered = expression.find("^+") != std::string::npos;
	}

	std::ostream& operator<<(std::ostream& os, const WickOperator& op)
	{
		os << "\\langle " << op.type << "_{ " << op.momentum << ", ";
		for (const auto& index : op.indizes) {
			os << index << " ";
		}
		os << "}";
		if (op.isDaggered) {
			os << "^\\dagger";
		}
		os << " \\rangle";
		return os;
	}
	std::ostream& operator<<(std::ostream& os, const std::vector<WickOperator>& ops)
	{
		for (const auto& op : ops) {
			os << op << " ";
		}
		return os;
	}
}
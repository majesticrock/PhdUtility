#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/utility/StringUtility.hpp>
#include <cassert>

namespace mrock::symbolic_operators {
	WickOperator::WickOperator(const OperatorType& _type, const bool _is_daggered, const Momentum& _momentum, const IndexWrapper& _indizes)
		: type(_type), is_daggered(_is_daggered), momentum(_momentum), indizes(_indizes) {}
	WickOperator::WickOperator(const OperatorType& _type, const bool _is_daggered, const Momentum& _momentum, const Index _index)
		: type(_type), is_daggered(_is_daggered), momentum(_momentum), indizes(_index) {}

	WickOperator::WickOperator(const std::string& expression)
	{
		// Syntax    type{Momentum_expression;index1,index2,...}(^+)

		this->type = string_to_wick.at(expression.substr(0U, expression.find('{')));
		std::vector<std::string> momentum_strings = mrock::utility::extract_elements(expression, '{', ';');
		std::vector<std::string> index_strings = mrock::utility::extract_elements(expression, ';', '}');

		assert(momentum_strings.size() == 1U);
		this->momentum = Momentum(momentum_strings.front());

		this->indizes.reserve(index_strings.size());
		for (const auto& arg : index_strings) {
			this->indizes.push_back(string_to_index.at(arg));
		}

		this->is_daggered = expression.find("^+") != std::string::npos;
	}

	std::ostream& operator<<(std::ostream& os, const WickOperator& op)
	{
		os << "\\langle " << op.type << "_{ " << op.momentum << ", ";
		for (const auto& index : op.indizes) {
			os << index << " ";
		}
		os << "}";
		if (op.is_daggered) {
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
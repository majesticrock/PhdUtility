#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/detail/string_helper.hpp>
#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/Operator.hpp>
#include <mrock/symbolic_operators/OperatorType.hpp>

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>
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
		std::vector<std::string> momentum_strings = extract_elements(expression, '{', ';');
		std::vector<std::string> index_strings = extract_elements(expression, ';', '}');

		assert(momentum_strings.size() == 1U);
		this->momentum = Momentum(momentum_strings.front());

		this->indizes.reserve(index_strings.size());
		for (const auto& arg : index_strings) {
			this->indizes.push_back(string_to_index.at(arg));
		}

		this->is_daggered = expression.find("^+") != std::string::npos;
	}

    std::vector<Operator> WickOperator::to_operator_expression() const
    {
		std::vector<Operator> result{
			Operator(this->momentum, this->indizes, false), 
			Operator(this->momentum, this->indizes, false)
		};

		switch(this->type) {
			case OperatorType::Eta:
				result[0].momentum.add_Q = true;
			case OperatorType::SC:
				result[0].momentum.flip_momentum();
				result[0].indizes.insert(result[0].indizes.begin(), Index::SpinDown);
				result[1].indizes.insert(result[1].indizes.begin(), Index::SpinUp);

				result[0].is_daggered = this->is_daggered;
				result[1].is_daggered = this->is_daggered;
				if(this->is_daggered) {
					std::swap(result[0], result[1]);
				}
				break;

			case OperatorType::CDW:
				result[0].momentum.add_Q = true;
			case OperatorType::Number:
				if(this->is_daggered) {
					std::swap(result[0], result[1]);
				}
				result[0].is_daggered = true;
				result[1].is_daggered = false;
				break;
				
			default:
				throw std::runtime_error("Operator type not handled!");
		}

		return result;
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
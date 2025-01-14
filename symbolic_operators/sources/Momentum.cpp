#include <mrock/symbolic_operators/Momentum.hpp>
#include <cctype>
#include <sstream>
#include <string>

namespace mrock::symbolic_operators {
	inline momentum_symbols::value_type identify_subexpression(const std::string& sub) {
		if (sub.front() == '+')
			return identify_subexpression(std::string(sub.begin() + 1, sub.end()));
		if (sub.front() == '-') {
			momentum_symbols::value_type ret = identify_subexpression(std::string(sub.begin() + 1, sub.end()));
			ret.factor *= -1;
			return ret;
		}
		if (!std::isdigit(sub.front()))
			return MomentumSymbol(1, sub.front());

		const auto it = std::find_if(sub.begin(), sub.end(), [](const char c) {
			return !std::isdigit(c);
			});

		return MomentumSymbol(std::stoi(std::string(sub.begin(), it)), sub.back());
	}

	Momentum::Momentum(const std::string& expression, bool Q/* = false*/) : add_Q(Q)
	{
		if(expression != "0") {
			size_t last = 0U;
			size_t current = expression.find_first_of("+-", expression.front() == '+' || expression.front() == '-' ? 1U : 0U);
			do {
				current = expression.find_first_of("+-", last + 1U);
				this->momentum_list.push_back(identify_subexpression(expression.substr(last, current - last)));
				last = current;
			} while (current != std::string::npos);
		}
	}

	void Momentum::sort()
	{
		for (size_t i = 0U; i < momentum_list.size(); ++i)
		{
			for (size_t j = i + 1U; j < momentum_list.size(); ++j)
			{
				// Comparing two chars is easy
				if (momentum_list[i].name > momentum_list[j].name) {
					std::swap(momentum_list[i], momentum_list[j]);
				}
			}
		}
		remove_zeros();
	}

	void Momentum::remove_contribution(const MomentumSymbol::name_type momentum) 
	{
		const int idx = this->isUsed(momentum);
		if (idx < 0) return;
		this->momentum_list.erase(this->momentum_list.begin() + idx);
	}

	Momentum& Momentum::operator+=(const Momentum& rhs)
	{
		this->add_Q = (rhs.add_Q != this->add_Q);
		bool foundOne = false;
		for (size_t i = 0U; i < rhs.momentum_list.size(); ++i)
		{
			foundOne = false;
			for (size_t j = 0U; j < this->momentum_list.size(); ++j)
			{
				if (rhs.momentum_list[i].name == this->momentum_list[j].name) {
					foundOne = true;
					this->momentum_list[j].factor += rhs.momentum_list[i].factor;
					if (this->momentum_list[j].factor == 0) {
						this->momentum_list.erase(this->momentum_list.begin() + j);
					}
					break;
				}
			}
			if (!foundOne) {
				this->momentum_list.push_back(rhs.momentum_list[i]);
			}
		}
		this->sort();
		return *this;
	}

	Momentum& Momentum::operator-=(const Momentum& rhs)
	{
		this->add_Q = (rhs.add_Q != this->add_Q);
		bool foundOne = false;
		for (size_t i = 0U; i < rhs.momentum_list.size(); ++i)
		{
			foundOne = false;
			for (size_t j = 0U; j < this->momentum_list.size(); ++j)
			{
				if (rhs.momentum_list[i].name == this->momentum_list[j].name) {
					foundOne = true;
					this->momentum_list[j].factor -= rhs.momentum_list[i].factor;
					if (this->momentum_list[j].factor == 0) {
						this->momentum_list.erase(this->momentum_list.begin() + j);
					}
					break;
				}
			}
			if (!foundOne) {
				this->momentum_list.push_back(MomentumSymbol(-rhs.momentum_list[i].factor, rhs.momentum_list[i].name));
			}
		}
		this->sort();
		return *this;
	}
	void Momentum::add_in_place(const Momentum& rhs)
	{
		(*this) += rhs;
	}

	void Momentum::replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith)
	{
		for (const auto& x : replaceWith.momentum_list) {
			if (x.name == replaceWhat) {
				throw std::invalid_argument("You are trying to replace a momentum with itself. This has undefined behaviour!");
			}
		}
		for (size_t i = 0U; i < momentum_list.size(); ++i) {
			if (momentum_list[i].name == replaceWhat) {
				auto buffer = replaceWith;
				buffer.multiply_by(momentum_list[i].factor);
				this->momentum_list.erase(momentum_list.begin() + i);

				(*this) += buffer;
			}
		}
	}

	void Momentum::remove_zeros()
	{
		for (auto it = momentum_list.begin(); it != momentum_list.end();) {
			if (it->factor == 0) {
				it = momentum_list.erase(it);
			}
			else {
				++it;
			}
		}
	}

	void Momentum::flip_single(const MomentumSymbol::name_type momentum)
	{
		for (auto& momentum_symbol : momentum_list) {
			if (momentum_symbol.name == momentum) {
				momentum_symbol.factor *= -1;
			}
		}
	}

	bool Momentum::operator==(const Momentum& rhs) const {
		if (this->add_Q != rhs.add_Q) return false;
		if (this->momentum_list.size() != rhs.momentum_list.size()) return false;
		bool foundOne = true;
		for (size_t i = 0U; i < this->momentum_list.size(); ++i)
		{
			foundOne = false;
			for (size_t j = 0U; j < rhs.momentum_list.size(); ++j)
			{
				if (this->momentum_list[i] == rhs.momentum_list[j]) {
					foundOne = true;
					break;
				}
			}
			if (!foundOne) return false;
		}
		return true;
	}

	std::string Momentum::to_string() const {
		std::ostringstream oss;
        oss << *this;
        return oss.str();
	}

	std::ostream& operator<<(std::ostream& os, const Momentum& momentum)
	{
		if (momentum.momentum_list.empty()) {
			if (momentum.add_Q) {
				os << "Q";
			}
			else {
				os << "0";
			}
			return os;
		}
		for (momentum_symbols::const_iterator it = momentum.momentum_list.begin(); it != momentum.momentum_list.end(); ++it)
		{
			if (it != momentum.momentum_list.begin() && it->factor > 0) {
				os << "+";
			}
			if (abs(it->factor) != 1) {
				os << it->factor;
			}
			else if (it->factor == -1) {
				os << "-";
			}
			os << it->name;
		}
		if (momentum.add_Q) {
			os << " + Q";
		}
		return os;
	}

	bool operator>(const Momentum& lhs, const Momentum& rhs)
	{
		if (lhs.momentum_list == rhs.momentum_list) return false;
		if (rhs.momentum_list.empty()) return true;
		if (lhs.momentum_list.empty()) return false;
		return lhs.momentum_list.front() > rhs.momentum_list.front();
	}
	bool operator<(const Momentum& lhs, const Momentum& rhs)
	{
		if (lhs.momentum_list == rhs.momentum_list) return false;
		if (lhs.momentum_list.empty()) return true;
		if (rhs.momentum_list.empty()) return false;
		return lhs.momentum_list.front() < rhs.momentum_list.front();
	}

	bool momentum_order(const Momentum& lhs, const Momentum& rhs)
	{
		if(rhs.momentum_list.empty()) {
			if(lhs.momentum_list.empty() && !lhs.add_Q && rhs.add_Q) return true;
			return false; 
		}
		if(lhs.momentum_list.empty()) return true;
		if(lhs.momentum_list[0].name < rhs.momentum_list[0].name) return true;
		if(lhs.momentum_list[0].name == rhs.momentum_list[0].name) {
			if(!lhs.add_Q && rhs.add_Q) return true;
		}
		return false;
	}
}
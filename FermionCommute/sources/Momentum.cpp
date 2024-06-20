#include "Momentum.hpp"
#include <cctype>

namespace SymbolicOperators {
	inline momentum_pairs::value_type identify_subexpression(const std::string& sub) {
		if (sub.front() == '+')
			return identify_subexpression(std::string(sub.begin() + 1, sub.end()));
		if (sub.front() == '-') {
			momentum_pairs::value_type ret = identify_subexpression(std::string(sub.begin() + 1, sub.end()));
			ret.first *= -1;
			return ret;
		}
		if (!std::isdigit(sub.front()))
			return std::make_pair(1, sub.front());

		const auto it = std::find_if(sub.begin(), sub.end(), [](const char c) {
			return !std::isdigit(c);
			});

		return std::make_pair(std::stoi(std::string(sub.begin(), it)), sub.back());
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
				if (momentum_list[i].second > momentum_list[j].second) {
					std::swap(momentum_list[i], momentum_list[j]);
				}
			}
		}
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
				if (rhs.momentum_list[i].second == this->momentum_list[j].second) {
					foundOne = true;
					this->momentum_list[j].first += rhs.momentum_list[i].first;
					if (this->momentum_list[j].first == 0) {
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
				if (rhs.momentum_list[i].second == this->momentum_list[j].second) {
					foundOne = true;
					this->momentum_list[j].first -= rhs.momentum_list[i].first;
					if (this->momentum_list[j].first == 0) {
						this->momentum_list.erase(this->momentum_list.begin() + j);
					}
					break;
				}
			}
			if (!foundOne) {
				this->momentum_list.push_back(std::make_pair(-rhs.momentum_list[i].first, rhs.momentum_list[i].second));
			}
		}
		this->sort();
		return *this;
	}
	void Momentum::addInPlace(const Momentum& rhs)
	{
		(*this) += rhs;
	}

	void Momentum::replaceOccurances(const char replaceWhat, const Momentum& replaceWith)
	{
		for (const auto& x : replaceWith.momentum_list) {
			if (x.second == replaceWhat) {
				throw std::invalid_argument("You are trying to replace a momentum with itself. This has undefined behaviour!");
			}
		}
		for (size_t i = 0U; i < momentum_list.size(); ++i) {
			if (momentum_list[i].second == replaceWhat) {
				auto buffer = replaceWith;
				buffer.multiplyMomentum(momentum_list[i].first);
				this->momentum_list.erase(momentum_list.begin() + i);

				(*this) += buffer;
			}
		}
	}

	void Momentum::remove_zeros()
	{
		for (auto it = momentum_list.begin(); it != momentum_list.end();) {
			if (it->first == 0) {
				it = momentum_list.erase(it);
			}
			else {
				++it;
			}
		}
	}

	void Momentum::flip_single(char momentum)
	{
		for (auto& momentum_pair : momentum_list) {
			if (momentum_pair.second == momentum) {
				momentum_pair.first *= -1;
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
		for (momentum_pairs::const_iterator it = momentum.momentum_list.begin(); it != momentum.momentum_list.end(); ++it)
		{
			if (it != momentum.momentum_list.begin() && it->first > 0) {
				os << "+";
			}
			if (abs(it->first) != 1) {
				os << it->first;
			}
			else if (it->first == -1) {
				os << "-";
			}
			os << it->second;
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
}
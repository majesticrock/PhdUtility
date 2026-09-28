#include <mrock/symbolic_operators/Exceptions.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace mrock::symbolic_operators {
// Private function used in string expression constructor
inline std::vector<MomentumSymbol>::value_type identify_subexpression(const std::string& sub) {
    if (sub.front() == '+')
        return identify_subexpression(std::string(sub.begin() + 1, sub.end()));
    if (sub.front() == '-') {
        std::vector<MomentumSymbol>::value_type ret = identify_subexpression(std::string(sub.begin() + 1, sub.end()));
        ret.factor *= -1;
        return ret;
    }
    if (!std::isdigit(sub.front()))
        return MomentumSymbol(1, sub.front());

    const auto it = std::find_if(sub.begin(), sub.end(), [](const char c) { return !std::isdigit(c); });

    return MomentumSymbol(std::stoi(std::string(sub.begin(), it)), sub.back());
}

void Momentum::sort() noexcept {
    std::sort(momentum_list.begin(), momentum_list.end());
    remove_zeros();
}

void Momentum::remove_contribution(const MomentumSymbol::name_type momentum) noexcept {
    const int idx = this->is_used_at(momentum);
    if (idx < 0)
        return;
    this->momentum_list.erase(this->momentum_list.begin() + idx);
}

void Momentum::add_in_place(const Momentum& rhs) {
    (*this) += rhs;
}

void Momentum::replace_occurances(const MomentumSymbol::name_type replaceWhat, const Momentum& replaceWith) {
    for (const auto& x : replaceWith.momentum_list) {
        if (x.name == replaceWhat) {
            throw momentum_replacement_error(static_cast<unsigned char>(replaceWhat), replaceWith.to_string());
        }
    }
    for (std::size_t i = 0U; i < momentum_list.size(); ++i) {
        if (momentum_list[i].name == replaceWhat) {
            auto buffer = replaceWith;
            buffer.multiply_by(momentum_list[i].factor);
            this->momentum_list.erase(momentum_list.begin() + i);

            (*this) += buffer;
        }
    }
    sort();
}

void Momentum::remove_zeros() noexcept {
    for (auto it = momentum_list.begin(); it != momentum_list.end();) {
        if (it->factor == 0) {
            it = momentum_list.erase(it);
        } else {
            ++it;
        }
    }
}

void Momentum::flip_single(const MomentumSymbol::name_type momentum) noexcept {
    for (auto& momentum_symbol : momentum_list) {
        if (momentum_symbol.name == momentum) {
            momentum_symbol.factor *= -1;
        }
    }
}

int Momentum::is_used_at(const MomentumSymbol::name_type value) const noexcept {
    for (std::size_t i = 0U; i < momentum_list.size(); ++i) {
        if (momentum_list[i].name == value)
            return i;
    }
    return -1;
}

std::string Momentum::to_string() const {
    std::ostringstream oss;
    oss << *this;
    return oss.str();
}

Momentum& Momentum::operator+=(const Momentum& rhs) {
    this->add_PI = (rhs.add_PI != this->add_PI);
    bool foundOne = false;
    for (std::size_t i = 0U; i < rhs.momentum_list.size(); ++i) {
        foundOne = false;
        for (std::size_t j = 0U; j < this->momentum_list.size(); ++j) {
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

Momentum& Momentum::operator-=(const Momentum& rhs) {
    this->add_PI = (rhs.add_PI != this->add_PI);
    bool foundOne = false;
    for (std::size_t i = 0U; i < rhs.momentum_list.size(); ++i) {
        foundOne = false;
        for (std::size_t j = 0U; j < this->momentum_list.size(); ++j) {
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

std::ostream& operator<<(std::ostream& os, const Momentum& momentum) {
    if (momentum.momentum_list.empty()) {
        if (momentum.add_PI) {
            os << _vector_wrap("\\Pi");
        } else {
            os << "0";
        }
        return os;
    }
    for (std::vector<MomentumSymbol>::const_iterator it = momentum.momentum_list.begin();
         it != momentum.momentum_list.end(); ++it) {
        if (it != momentum.momentum_list.begin() && it->factor > 0) {
            os << "+";
        }
        if (std::abs(it->factor) != 1) {
            os << it->factor;
        } else if (it->factor == -1) {
            os << "-";
        }
        os << _vector_wrap(it->name);
    }
    if (momentum.add_PI) {
        os << " + " + _vector_wrap("\\Pi");
    }
    return os;
}

Momentum::Momentum(const char value, int plus_minus /* = 1 */, bool add_PI_ /* = false */)
    : momentum_list(1, MomentumSymbol(plus_minus, value)), add_PI(add_PI_) 
{
    sort();
}

Momentum::Momentum(const MomentumSymbol::name_type value, int plus_minus /* = 1 */, bool add_PI_ /* = false */)
    : momentum_list(1, {plus_minus, value}), add_PI(add_PI_) 
{
    sort();
}

Momentum::Momentum(const std::vector<MomentumSymbol>& _momenta, bool add_PI_ /* = false */)
    : momentum_list(_momenta), add_PI(add_PI_) 
{
    sort();
}

Momentum::Momentum(MomentumSymbol const& momentum_symbol, bool add_PI_ /* = false */)
    : momentum_list{momentum_symbol}, add_PI(add_PI_) 
{
    sort();
}

Momentum::Momentum(const std::string& expression, bool add_PI_ /* = false*/) : add_PI(add_PI_) {
    if (expression != "0") {
        std::size_t last = 0U;
        std::size_t current =
            expression.find_first_of("+-", expression.front() == '+' || expression.front() == '-' ? 1U : 0U);
        do {
            current = expression.find_first_of("+-", last + 1U);
            this->momentum_list.push_back(identify_subexpression(expression.substr(last, current - last)));
            last = current;
        } while (current != std::string::npos);
    }
    sort();
}

std::strong_ordering Momentum::operator<=>(const Momentum& other) const noexcept {
    if (auto cmp = momentum_list.size() <=> other.momentum_list.size(); cmp != 0)
        return cmp;
    if (auto cmp = momentum_list <=> other.momentum_list; cmp != 0)
        return cmp;
    
    return add_PI <=> other.add_PI;
}

}  // namespace mrock::symbolic_operators
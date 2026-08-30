#include <mrock/symbolic_operators/Coefficient.hpp>
#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/MomentumList.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>
#include <mrock/symbolic_operators/detail/string_helper.hpp>

#include <map>
#include <ostream>
#include <utility>

namespace mrock::symbolic_operators {
void Coefficient::invert_momentum(const MomentumSymbol::name_type what) {
    for (auto& mom : momenta) {
        if (this->inversion_symmetry && mom.size() == 1U) {
            // If the coefficient is symmetric under inversion, i.e., c(k) = c(-k),
            // we want to make sure that the momentum is always positive
            if (mom.front().factor > 0)
                continue;
        }
        mom.flip_single(what);
    }
}

void Coefficient::use_symmetric_interaction_exchange() {
    if (this->is_symmetrized_interaction) {
        assert(this->momenta.size() == 3U);
        std::swap(this->momenta[0], this->momenta[1]);
        this->momenta.back().flip_momentum();
    }
}

void Coefficient::use_symmetric_interaction_inversion() {
    if (this->is_symmetrized_interaction) {
        assert(this->momenta.size() == 3U);
        for (auto& mom : momenta) {
            mom.flip_momentum();
        }
    }
}

void Coefficient::use_inversion_symmetry() {
    if (!this->inversion_symmetry)
        return;
    for (auto& momentum : momenta) {
        momentum.flip_momentum();
    }
}

void Coefficient::remove_momentum_contribution(const MomentumSymbol::name_type value) {
    for (auto& mom : momenta) {
        mom.remove_contribution(value);
    }
}

void Coefficient::use_custom_symmetry() {
    if (custom_symmetry.has_value()) {
        custom_symmetry.value()(*this);
    }
}

Coefficient Coefficient::parse_string(const std::string& expression,
                                      bool _Q_changes_sign /* = false */,
                                      bool _inversion_symmetry /* = true */) {
    // Syntax:   name{Momentum_expression1,Momentum_expression1;index1,index2,...}
    Coefficient ret;
    ret.name = expression.substr(0U, find_skip_escaped(expression, '{'));
    remove_escape_characters(ret.name);
    std::vector<std::string> momentum_strings = extract_elements(expression, '{', ';');
    std::vector<std::string> index_strings = extract_elements(expression, ';', '}');

    ret.momenta.reserve(momentum_strings.size());
    for (const auto& arg : momentum_strings) {
        ret.momenta.push_back(Momentum(arg));
    }
    ret.indices.reserve(index_strings.size());
    for (const auto& arg : index_strings) {
        ret.indices.push_back(string_to_index.at(arg));
    }

    ret.Q_changes_sign = _Q_changes_sign;
    ret.inversion_symmetry = _inversion_symmetry;
    return ret;
}

Coefficient Coefficient::parse_interaction_string(const std::string& expression) {
    Coefficient ret = parse_string(expression, false, false);
    ret.is_symmetrized_interaction = true;
    return ret;
}

Coefficient::Coefficient(const std::string& _name)
    : name(_name), momenta(), indices(), Q_changes_sign(false), is_daggered(false) {}

Coefficient::Coefficient(const std::string& _name,
                         const Momentum& _momentum,
                         const IndexWrapper& _indices,
                         bool _Q_changes_sign,
                         bool _inversion_symmetry,
                         bool _is_daggered)
    : name(_name),
      momenta(_momentum),
      indices(_indices),
      inversion_symmetry{_inversion_symmetry},
      Q_changes_sign(_Q_changes_sign),
      is_daggered(_is_daggered) {}

Coefficient::Coefficient(const std::string& _name,
                         const Momentum& _momentum,
                         bool _Q_changes_sign,
                         bool _inversion_symmetry,
                         bool _is_daggered)
    : name(_name),
      momenta(_momentum),
      indices(),
      inversion_symmetry{_inversion_symmetry},
      Q_changes_sign(_Q_changes_sign),
      is_daggered(_is_daggered) {}

Coefficient::Coefficient(const std::string& _name,
                         const MomentumList& _momenta,
                         const IndexWrapper& _indices,
                         bool _Q_changes_sign,
                         bool _inversion_symmetry,
                         bool _is_daggered)
    : name(_name),
      momenta(_momenta),
      indices(_indices),
      inversion_symmetry{_inversion_symmetry},
      Q_changes_sign(_Q_changes_sign),
      is_daggered(_is_daggered) {}

Coefficient Coefficient::RealInversionSymmetric(
    const std::string& name,
    const MomentumList& momenta,
    const std::optional<std::function<void(Coefficient&)>>& custom_symmetry /* = std::nullopt */) {
    Coefficient ret(name, momenta, {}, true, false);
    ret.custom_symmetry = custom_symmetry;
    return ret;
}
Coefficient Coefficient::RealInteraction(
    const std::string& name,
    const MomentumList& momenta,
    const std::optional<std::function<void(Coefficient&)>>& custom_symmetry /* = std::nullopt */) {
    assert(momenta.size() == 3U);
    Coefficient ret(name, momenta, {}, false, false);
    ret.is_symmetrized_interaction = true;
    ret.custom_symmetry = custom_symmetry;
    return ret;
}

Coefficient Coefficient::RealInteraction(const std::string& name,
                                         const MomentumList& momenta,
                                         const IndexWrapper& _indices,
                                         const std::optional<std::function<void(Coefficient&)>>& custom_symmetry) {
    assert(momenta.size() == 3U);
    Coefficient ret(name, momenta, _indices, false, false);
    ret.is_symmetrized_interaction = true;
    ret.custom_symmetry = custom_symmetry;
    return ret;
}

Coefficient Coefficient::HoneyComb(
    const std::string& name,
    const Momentum& momentum,
    bool daggered,
    bool is_real /* = true */,
    const std::optional<std::function<void(Coefficient&)>>& custom_symmetry /* = std::nullopt */) {
    Coefficient ret(name, momentum, {}, false, false, daggered);
    ret.custom_symmetry = custom_symmetry;
    ret.is_real = is_real;
    return ret;
}

Coefficient Coefficient::Constant(const std::string& name,
                                  const IndexWrapper& indices /* = IndexWrapper{} */,
                                  bool is_real /* = false */,
                                  bool is_daggered /* = true */) {
    Coefficient ret(name, MomentumList{}, indices, false, true, is_daggered);
    ret.is_real = is_real;
    return ret;
}

bool Coefficient::uses_index(const Index index) const noexcept {
    for (const auto& idx : indices) {
        if (idx == index)
            return true;
    }
    return false;
}

bool Coefficient::depends_on_momentum() const noexcept {
    if (this->momenta.empty())
        return false;
    return std::any_of(this->momenta.begin(), this->momenta.end(),
                       [](const Momentum& momentum) { return !momentum.momentum_list.empty(); });
}

bool Coefficient::depends_on(const MomentumSymbol::name_type momentum) const noexcept {
    if (this->momenta.empty())
        return false;
    return std::any_of(this->momenta.begin(), this->momenta.end(),
                       [momentum](const Momentum& mom) { return mom.uses(momentum); });
}

// This function determines whether the coefficient depends on something like k-l
// Currently, this only makes sense if the coefficient does not depend on
bool Coefficient::depends_on_two_momenta() const noexcept {
    assert(momenta.size() == 1U);
    return this->momenta.front().momentum_list.size() == 2U;
}

Coefficient& Coefficient::hermitian_conjugate_inplace() {
    if (is_real) {
        is_daggered = false;
        return *this;
    }
    is_daggered = !is_daggered;
    return *this;
}

Coefficient Coefficient::hermitian_conjugate() const {
    Coefficient copy(*this);
    copy.hermitian_conjugate_inplace();
    return copy;
}

std::ostream& operator<<(std::ostream& os, const Coefficient& coeff) {
    os << coeff.name;
    if (!coeff.indices.empty()) {
        os << "_{ " << coeff.indices << "}";
    }
    if (coeff.is_daggered) {
        os << "^*";
    }
    os << coeff.momenta << " ";
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Coefficient>& coeffs) {
    for (auto& coeff : coeffs) {
        os << coeff << " ";
    }
    return os;
}

bool operator==(const Coefficient& lhs, const Coefficient& rhs) {
    if (lhs.name != rhs.name)
        return false;
    if (lhs.momenta != rhs.momenta)
        return false;
    if (lhs.is_daggered != rhs.is_daggered)
        return false;
    return (lhs.indices == rhs.indices);
}

}  // namespace mrock::symbolic_operators
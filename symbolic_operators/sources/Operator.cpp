#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/Operator.hpp>

#include <ostream>

namespace mrock::symbolic_operators {
std::ostream& operator<<(std::ostream& os, const Operator& op) {
    if (op.is_fermion) {
        os << "\\hat{c}";
    } else {
        os << "\\hat{b}";
    }
    os << "_{ " << op.momentum << ", " << op.indices << "}";
    if (op.is_daggered) {
        os << "^\\dagger ";
    } else {
        os << "^{\\phantom{\\dagger}}";  // For alignment purposes
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Operator>& ops) {
    for (const auto& op : ops) {
        os << op;
    }
    return os;
}

std::strong_ordering Operator::operator<=>(const Operator& other) const {
    if (auto cmp = is_fermion <=> other.is_fermion; cmp != 0)
        return cmp;
    if (auto cmp = is_daggered <=> other.is_daggered; cmp != 0)
        return cmp;
    if (auto cmp = indices <=> other.indices; cmp != 0)
        return cmp;
    return momentum <=> other.momentum;
}

Operator::Operator(const Momentum& _momentum, const IndexWrapper _indices, bool _is_daggered, bool _is_fermion)
    : momentum(_momentum), indices(_indices), is_daggered(_is_daggered), is_fermion(_is_fermion) {}

Operator::Operator(const std::vector<MomentumSymbol>& _momentum,
                   const IndexWrapper _indices,
                   bool _is_daggered,
                   bool _is_fermion)
    : momentum(_momentum), indices(_indices), is_daggered(_is_daggered), is_fermion(_is_fermion) {}

Operator::Operator(const MomentumSymbol::name_type _momentum,
                   bool add_PI,
                   const IndexWrapper _indices,
                   bool _is_daggered,
                   bool _is_fermion)
    : momentum(_momentum, add_PI), indices(_indices), is_daggered(_is_daggered), is_fermion(_is_fermion) {}

Operator::Operator(const MomentumSymbol::name_type _momentum,
                   int sign,
                   bool add_PI,
                   const IndexWrapper _indices,
                   bool _is_daggered,
                   bool _is_fermion)
    : momentum(_momentum, sign, add_PI), indices(_indices), is_daggered(_is_daggered), is_fermion(_is_fermion) {}
}  // namespace mrock::symbolic_operators
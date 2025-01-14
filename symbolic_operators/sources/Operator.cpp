#include <mrock/symbolic_operators/Operator.hpp>

namespace mrock::symbolic_operators {
	std::ostream& operator<<(std::ostream& os, const Operator& op)
	{
		if (op.is_fermion) {
			os << "c";
		} else {
			os << "b";
		}
		os << "_{ " << op.momentum << ", " << op.indizes << "}";
		if (op.is_daggered) {
			os << "^\\dagger ";
		}
		return os;
	}

	std::ostream& operator<<(std::ostream& os, const std::vector<Operator>& ops)
	{
		for (const auto& op : ops) {
			os << op;
		}
		return os;
	}

	Operator::Operator(const Momentum& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion)
		: momentum(_momentum), indizes(_indizes), is_daggered(_is_daggered), is_fermion(_is_fermion) {}

	Operator::Operator(const momentum_symbols& _momentum, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion)
		: momentum(_momentum), indizes(_indizes), is_daggered(_is_daggered), is_fermion(_is_fermion) {}

	Operator::Operator(const MomentumSymbol::name_type _momentum, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion)
		: momentum(_momentum, add_Q), indizes(_indizes), is_daggered(_is_daggered), is_fermion(_is_fermion) {}

	Operator::Operator(const MomentumSymbol::name_type _momentum, int sign, bool add_Q, const IndexWrapper _indizes, bool _is_daggered, bool _is_fermion)
		: momentum(_momentum, sign, add_Q), indizes(_indizes), is_daggered(_is_daggered), is_fermion(_is_fermion) {}
}
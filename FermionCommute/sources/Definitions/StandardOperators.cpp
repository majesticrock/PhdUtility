#include "StandardOperators.hpp"

namespace SymbolicOperators {
	const Momentum StandardOperators::base_k{ momentum_pairs{ {1, 'k'} }, false };
	const Momentum StandardOperators::base_x{ momentum_pairs{ {1, 'x'} }, false };

	const Operator StandardOperators::c_k{ base_k, SpinUp, false };
	const Operator StandardOperators::c_minus_k{ -base_k, SpinDown, false };
	const Operator StandardOperators::c_k_dagger{ base_k, SpinUp, true };
	const Operator StandardOperators::c_minus_k_dagger{ -base_k, SpinDown, true };
}
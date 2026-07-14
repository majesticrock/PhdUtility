#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

namespace mrock::symbolic_operators {
void SpinSymmetry::apply_to(WickTerm& term) const {
    for (auto& op : term.operators) {
        for (auto& idx : op.indizes) {
            if (idx == Index::SpinDown)
                idx = Index::SpinUp;
        }
    }
}

void InversionSymmetry::apply_to(WickTerm& term) const {
    for (auto& op : term.operators) {
        if (!op.momentum.momentum_list.empty() && op.momentum.momentum_list[0].factor < 0) {
            op.momentum.flip_momentum();
        }
    }
}
}  // namespace mrock::symbolic_operators
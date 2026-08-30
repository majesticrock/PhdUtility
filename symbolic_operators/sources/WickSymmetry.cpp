#include <mrock/symbolic_operators/IndexWrapper.hpp>
#include <mrock/symbolic_operators/Momentum.hpp>
#include <mrock/symbolic_operators/MomentumSymbol.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

namespace mrock::symbolic_operators {
void SpinSymmetry::apply_to(std::vector<WickOperator>& operator_vector) const {
    for (auto& op : operator_vector) {
        for (auto& idx : op.indices) {
            if (idx == Index::SpinDown)
                idx = Index::SpinUp;
        }
    }
}

void InversionSymmetry::apply_to(std::vector<WickOperator>& operator_vector) const {
    for (auto& op : operator_vector) {
        if (!op.momentum.momentum_list.empty() && op.momentum.momentum_list[0].factor < 0) {
            op.momentum.flip_momentum();
        }
    }
}
}  // namespace mrock::symbolic_operators
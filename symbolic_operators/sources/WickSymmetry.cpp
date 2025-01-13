#include <mrock/symbolic_operators/WickSymmetry.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

namespace mrock::symbolic_operators{
    void SpinSymmetry::apply_to(WickTerm& term) const
    {
        for (auto& op : term.operators) {
            for (auto& idx : op.indizes) {
                if(idx == Index::SpinDown)
                    idx = Index::SpinUp;
            }
        }
    }

    void InversionSymmetry::apply_to(WickTerm& term) const
    {
        for (auto& op : term.operators) {
            if (!op.momentum.momentum_list.empty() && op.momentum.momentum_list[0].first < 0) {
                op.momentum.flip_momentum();
            }
        }
    }
}
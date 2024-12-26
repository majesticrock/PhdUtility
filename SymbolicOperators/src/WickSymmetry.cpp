#include <mrock/SymbolicOperators/WickSymmetry.hpp>
#include <mrock/SymbolicOperators/WickTerm.hpp>

namespace mrock::SymbolicOperators{
    void SpinSymmetry::apply_to(WickTerm& term) const
    {
        for (auto& op : term.operators) {
            for (auto& idx : op.indizes) {
                if(idx == Index::SpinDown)
                    idx = Index::SpinUp;
            }
        }
    }

    void TranslationalSymmetry::apply_to(WickTerm& term) const
    {
        for (auto& op : term.operators) {
            if (!op.momentum.momentum_list.empty() && op.momentum.momentum_list[0].first < 0) {
                op.momentum.flipMomentum();
            }
        }
    }
}
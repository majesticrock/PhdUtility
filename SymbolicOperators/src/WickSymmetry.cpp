#include <SymbolicOperators/WickSymmetry.hpp>
#include <SymbolicOperators/WickTerm.hpp>

namespace SymbolicOperators{
    void SpinSymmetry::apply_to(WickTerm& term) const
    {
        for (auto& op : term.operators) {
            for(auto& idx : op.indizes){
                if(idx == SpinDown)
                    idx = SpinUp;
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
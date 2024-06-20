#pragma once
#include "OperatorType.hpp"
#include "WickTerm.hpp"
#include <type_traits>

namespace SymbolicOperators
{
    struct WickSymmetry {
        virtual void apply_to(WickTerm& term) const = 0;
    };

    // Expectation values for spin up and down are the same
    struct SpinSymmetry : public WickSymmetry {
        void apply_to(WickTerm& term) const override;
    };
    // Expectation values for k and -k are the same
    struct TranslationalSymmetry : public WickSymmetry {
        void apply_to(WickTerm& term) const override;
    };
    // <operator^+> = <operator>
    template<OperatorType... operators>
    struct PhaseSymmetry : public WickSymmetry {
        void apply_to(WickTerm& term) const override
        {
            for(auto& op : term.operators) {
                auto equals_op = [&op](OperatorType comp){
                    return op.type == comp;
                    };
                if ( (equals_op(operators) || ...) )
                {
                    op.isDaggered = false;
                }
            }
        };
    };
}

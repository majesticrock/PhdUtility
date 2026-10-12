#include <mrock/symbolic_operators/WickOperator.hpp>

#include <cassert>
#include <string>

using namespace mrock::symbolic_operators;

int main() {
    const Momentum k('k');
    const WickOperator original{OperatorType::SC, true, k, IndexWrapper{Index::SpinUp}};

    const Coefficient coeff = original.to_coefficient("Delta");
    assert(coeff.name == "Delta_{f}");
    assert(coeff.momenta.size() == 1U);
    assert(coeff.momenta.front() == k);
    assert(coeff.indices == IndexWrapper{Index::SpinUp});
    assert(coeff.is_daggered);

    const WickOperator roundtrip = WickOperator::from_coefficient(coeff, "Delta");
    assert(roundtrip.type == OperatorType::SC);
    assert(roundtrip.is_daggered == original.is_daggered);
    assert(roundtrip.momentum == original.momentum);
    assert(roundtrip.indices == original.indices);

    const Coefficient unmatched{"Other", Momentum('p')};
    const WickOperator invalid = WickOperator::from_coefficient(unmatched, "Delta");
    assert(invalid.type == OperatorType::Undefined);

    return 0;
}

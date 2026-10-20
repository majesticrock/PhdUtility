#include <mrock/symbolic_operators/WickOperator.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

using namespace mrock::symbolic_operators;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_throws_for_momentum_count(const Coefficient& coefficient) {
    try {
        (void)WickOperator::from_coefficient(coefficient, "Delta");
    } catch (const std::runtime_error&) {
        return;
    }
    throw std::runtime_error("from_coefficient accepted a coefficient without exactly one momentum");
}

void test_round_trip_for_each_type() {
    const std::string base = "Delta";
    const Momentum momentum('k');
    const IndexWrapper indices{Index::SpinUp, Index::SpinDown};

    for (const auto& [type, type_name] : wick_to_string) {
        for (const bool daggered : {false, true}) {
            const WickOperator original{type, daggered, momentum, indices};
            const Coefficient coefficient = original.to_coefficient(base);

            require(coefficient.name == base + "_{(" + type_name + ")}", "to_coefficient generated an unexpected name");
            require(coefficient.momenta.size() == 1U && coefficient.momenta.front() == momentum,
                    "to_coefficient did not preserve momentum");
            require(coefficient.indices == indices, "to_coefficient did not preserve indices");
            require(coefficient.is_daggered == daggered, "to_coefficient did not preserve dagger state");
            require(!coefficient.custom_symmetry.has_value(), "to_coefficient retained a custom symmetry");
            require(!coefficient.inversion_symmetry && !coefficient.is_real && !coefficient.is_symmetrized_interaction,
                    "to_coefficient did not clear coefficient symmetries");

            const WickOperator round_trip = WickOperator::from_coefficient(coefficient, base);
            require(round_trip.type == type, "from_coefficient did not restore operator type");
            require(round_trip.is_daggered == daggered, "from_coefficient did not restore dagger state");
            require(round_trip.momentum == momentum, "from_coefficient did not restore momentum");
            require(round_trip.indices == indices, "from_coefficient did not restore indices");

            Coefficient legacy{base + "{(" + type_name + ")}", momentum, indices};
            legacy.is_daggered = daggered;
            const WickOperator legacy_round_trip = WickOperator::from_coefficient(legacy, base);
            require(legacy_round_trip == original, "from_coefficient rejected the legacy coefficient name");
        }
    }
}

void test_unmatched_names() {
    const Coefficient wrong_base{"Other_{(n)}", Momentum('p')};
    const Coefficient wrong_suffix{"Delta_{(unknown)}", Momentum('p')};
    const Coefficient partial_base{"DeltaExtra_{(n)}", Momentum('p')};

    for (const Coefficient* coefficient : {&wrong_base, &wrong_suffix, &partial_base}) {
        const WickOperator unmatched = WickOperator::from_coefficient(*coefficient, "Delta");
        require(unmatched.type == OperatorType::Undefined, "from_coefficient matched an unrelated name");
        require(unmatched.momentum == coefficient->momenta.front(), "unmatched conversion did not preserve momentum");
    }
}

void test_invalid_momentum_counts() {
    require_throws_for_momentum_count(Coefficient::Constant("Delta_{(n)}"));
    require_throws_for_momentum_count(Coefficient{"Delta_{(n)}", MomentumList{Momentum('k'), Momentum('p')}});
}

}  // namespace

int main() {
    try {
        test_round_trip_for_each_type();
        test_unmatched_names();
        test_invalid_momentum_counts();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}

#include <mrock/symbolic_operators/SumContainer.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>

#include <algorithm>

namespace mrock::symbolic_operators {
SumContainer& SumContainer::append(const SumContainer& other) {
    append_vector(this->momenta, other.momenta);
    append_vector(this->spins, other.spins);
    return *this;
}

SumContainer& SumContainer::append(const MomentumSum& other) {
    append_vector(this->momenta, other);
    return *this;
}

SumContainer& SumContainer::append(const IndexSum& other) {
    append_vector(this->spins, other);
    return *this;
}

SumContainer& SumContainer::prepend(const SumContainer& other) {
    prepend_vector(this->momenta, other.momenta);
    prepend_vector(this->spins, other.spins);
    return *this;
}

SumContainer& SumContainer::prepend(const MomentumSum& other) {
    prepend_vector(this->momenta, other);
    return *this;
}

SumContainer& SumContainer::prepend(const IndexSum& other) {
    prepend_vector(this->spins, other);
    return *this;
}

void SumContainer::sort() {
    std::sort(momenta.begin(), momenta.end());
    std::sort(spins.begin(), spins.end());
}

}  // namespace mrock::symbolic_operators

#include <mrock/symbolic_operators/SumContainer.hpp>
#include <mrock/symbolic_operators/detail/container_helper.hpp>

namespace mrock::symbolic_operators {
    SumContainer &mrock::symbolic_operators::SumContainer::append(const SumContainer &other)
    {
        append_vector(this->momenta, other.momenta);
    	append_vector(this->spins, other.spins);
    	return *this;
    }

    SumContainer &mrock::symbolic_operators::SumContainer::append(const MomentumSum &other)
    {
        append_vector(this->momenta, other);
		return *this;
    }

    SumContainer &mrock::symbolic_operators::SumContainer::append(const IndexSum &other)
    {
        append_vector(this->spins, other);
		return *this;
    }
}

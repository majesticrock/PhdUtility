#include <mrock/symbolic_operators/SumContainer.hpp>

namespace mrock::symbolic_operators {
    SumContainer &mrock::symbolic_operators::SumContainer::append(const SumContainer &other)
    {
        mrock::utility::append_vector(this->momenta, other.momenta);
    	mrock::utility::append_vector(this->spins, other.spins);
    	return *this;
    }

    SumContainer &mrock::symbolic_operators::SumContainer::append(const MomentumSum &other)
    {
        mrock::utility::append_vector(this->momenta, other);
		return *this;
    }

    SumContainer &mrock::symbolic_operators::SumContainer::append(const IndexSum &other)
    {
        mrock::utility::append_vector(this->spins, other);
		return *this;
    }
}

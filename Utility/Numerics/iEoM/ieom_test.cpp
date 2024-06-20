#include "_internal_functions.hpp"
#include <iostream>

//using RealType = double;
template<class Derived, class RealType>
struct ResTest {
public:
	void res_func() {
		_derived->createStartingStates();
		_derived->fillMatrices();

		std::cout << "res_func()" << std::endl;

		Eigen::VectorXd c = Eigen::VectorXd::Zero(10);
		_internal.template applyMatrixOperation<Utility::Numerics::iEoM::IEOM_INVERSE>(c);
	}

	ResTest(Derived* derived_ptr, RealType const& sqrt_precision)
		: _internal(sqrt_precision), _derived(derived_ptr) { };

public:
	Utility::Numerics::iEoM::ieom_internal<RealType> _internal;
	Derived* _derived;
};

struct XPTest : private ResTest<XPTest, double> {
	using parent = ResTest<XPTest, double>;
	//friend struct ResTest<XPTest>;

	void fill_M() {
		std::cout << "fill_M()" << std::endl;
	}
	void createStartingStates() {
		std::cout << "createStartingStates()" << std::endl;
	}
	void fillMatrices() {
		std::cout << "fillMatrices()" << std::endl;
	};

	void computeCollectiveModes() {
		std::cout << "computeCollectiveModes()" << std::endl;
		parent::res_func();
	}

	XPTest() : parent(this, 1e-6) {};
};

int main(int argc, char** argv) {
	XPTest test;
	test.computeCollectiveModes();

	return 0;
}
# Using the iEoM as described in https://arxiv.org/abs/2403.18701 to obtain excitation spectra

## XPResolvent.hpp and GeneralResolvent.hpp
Both work in essentially the same way.
A user-defined class should inherit from either of the algorithm-providing classes in the following way

```
class MyClass : public mrock::utility::Numerics::iEoM::XPResolvent<MyClass, double>
```

The user is required to define the following methods:

### fill_M()
Needs to fill the dynamical matrix *K_plus / K_minus* for ```XPResolvent``` and *M* for ```GeneralResolvent```

### fillMatrices()
Needs to fill both the dynamical matrix (see above) and norm matrix *N*

### createStartingStates()
Needs to create the starting states for the Lanczos algorithm.
These vectors govern which kind of Green's function will be computed.
For ```XPResolvent``` the starting states are each an ```std::array<Vector, 2>```, 
where the first (second) entry corresponds to *K_plus (K_minus)*.
For ```GeneralResolvent``` the starting states are plain vectors.

Naturally, the classes need access to these methods either by being a friend or the methods being publically visible.

Then ```compute_collective_modes``` and ```dynamic_matrix_is_negative``` provide the algorithms to either compute the excitations of the system or to simply determine whether *M* is non-negative.
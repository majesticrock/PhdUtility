# mrock/utility

This sublibrary contains various functionality that is required by the other applications in the larger superproject, such as `Hubbard`, `LatticeCUT`, etc.
The `iEoM` and `symbolic_operators` are designed to be independent of each other and of this sublibrary so that they can work as standalone libraries.
Thus, if you are just interested in using those, you can skip this subproject entirely.

That being said, you can of course reuse the functionalities of this subproject where needed.

Here, I will provide a list of functionalities, but will not go into as much detail as for the other two libraries.
The tests provided in the `tests` subdir serve as examples of how to use various scripts.

## Prerequisites

These prerequisites do not apply to every script, but if you want to use everything, you will need them.

- Eigen https://eigen.tuxfamily.org/index.php?title=Main_Page or https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost https://www.boost.org/
- nlohmann/json.hpp  https://json.nlohmann.me/ or https://github.com/nlohmann/json

## Header overview

### better_to_string.hpp
Provides fixed-precision to_string functionality.

### BinaryIO.hpp
Helps with writing and reading binary files.

### ComplexNumberIterators.hpp
Allows iterating just the real or just the imaginary part of a container of `std::complex` numbers.

### constexpr_power.hpp
Evaluates $n^m$ at compile time if $n$ and $m$ are compile-time constants.

### function_time.hpp
Measures the time it takes to run a function.

### info_to_json.hpp
The mrock CMake functionality creates a metadata header file.
The `info_to_json.hpp` header turns the `char` arrays in that metadata file into a `nlohmann::json` object.

### InputFileReader.hpp
Provides functionality to read config files; this is extensively used in `Hubbard`, `LatticeCUT`, and `ContinuumSystem`.

### is_complex.hpp
Provides a type-traits-like struct that asks whether or not a given type is a version of `std::complex`.

### OutputConvenience.hpp and OutputWriter.hpp
Provides helper functionality for writing data to files.
Unless the preprocessor flag `_NO_BOOST` is defined, it will use boost's serialization features and save the data as zip-compressed objects.

### progress_bar.hpp
Provides a rudimentary progress bar. An example is provided in `apps/progress_bar.cpp`.

### ThrowException.hpp
Throws an exception only if `NDEBUG` is not defined.

### UnderlyingRealType.hpp
Provides functionality that returns the template argument of `std::complex<Real>` or `Real` if `Real` is not an instance of `std::complex`.


## Numerics

### Integration/AdaptiveTrapezoidalRule.hpp
Provides an implementation for an adaptive-step-size trapezoidal rule.

### Integration/CauchyPrincipalValue.hpp and GeneralizedPrincipalValue.hpp
Provides functionality to compute principal-value integrals.

### Roots/Bisection.hpp
Implements the bisection method for finding roots of a one-dimensional function.

### Roots/BroydensMethod.hpp
Implements Broyden's method to find a root of a multivariate function.

### ErrorFunctors.hpp
Used to gauge the errors of numerical methods.

### hypergeometric_2F1.hpp
Computes the hypergeometric 2F1 function, see https://en.wikipedia.org/wiki/Hypergeometric_function or https://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/
This is used in `LatticeCUT` to evaluate certain densities of states.

### Interpolation.hpp
Provides interpolation functionality.

### is_integer.hpp
Asks whether a floating-point number is within a given tolerance of an integer.



## Selfconsistency

### Prerequisites

Some model class which defines

```
void iteration_step(const Eigen::Vector& input, const Eigen::Vector& error)
where error = input - new_delta

std::string info() const
```

and some attribute class that defines

```
bool converged
begin()
end()
size()
operator[]
```

These classes are used by `Hubbard`, `LatticeCUT` and `ContinuumSystem`.

### BroydenSolver.hpp
Uses `Roots/BroydensMethod.hpp` to solve a self-consistency problem.

### IterativeSolver.hpp
Solves a self-consistency problem iteratively.
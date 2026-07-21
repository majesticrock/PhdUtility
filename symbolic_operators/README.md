# symbolic_operators

`symbolic_operators` is a lightweight C++20 library for symbolic manipulation of canonical creation and annihilation operators. 

It supports both fermionic and bosonic operators and can reduce higher-order expectation values to bilinear products using Wick’s theorem. 

## Requirements

Required dependencies are a C++20 compiler and a functioning build environment. 

CMake 3.30 or newer is recommended for easier installation and running tests. 

Boost is optional and can be used for serialization, saving, and reading results. 

## Basic Usage

User-facing functionality is defined in the namespace `mrock::symbolic_operators`. 

Typical workflows create symbolic expressions with `Term`, manipulate them with functions such as `commutator`, `normal_order`, and `clean_up`, and then apply `wicks_theorem` to obtain `WickTerm` objects. 

```cpp
std::vector<Term> commutator_result = commutator(H, right);
clean_up(commutator_result);

WickTermCollector wicks;
wicks_theorem(commutator_result, templates, wicks);
clean_wicks(wicks);
```

## Main Classes

- `Term`: represents symbolic strings of creation and annihilation operators. 
- `Operator`: represents a single fermionic or bosonic creation or annihilation operator. 
- `Coefficient`: represents symbolic coefficients such as hopping amplitudes or interactions. 
- `Momentum`, `MomentumSymbol`, and `MomentumList`: provide symbolic momentum algebra. 
- `SumContainer` and `SymbolicSum`: store symbolic momentum and index sums. 
- `WickOperatorTemplate`: defines contraction patterns recognized by Wick’s theorem. 
- `WickTerm` and `WickTermCollector`: store Wick-expanded expressions and collections of such terms. 
- `WickSymmetry`: base class for simplifying Wick-expanded expressions using symmetries such as spin, inversion, or phase symmetry. 

## Example Workflow

A compact BCS-style example is available in `tests/bcs.cpp`. 

The example builds a Hamiltonian, computes symbolic expressions such as commutators, applies Wick’s theorem, simplifies the result with symmetries, and evaluates the expressions numerically. 

## More Information

For the complete API overview, class attributes, Wick-contraction setup, symmetry handling, and the full example workflow, see `user_guide_symop`. 
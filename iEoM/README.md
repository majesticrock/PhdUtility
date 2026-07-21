# iEoM

`iEoM` is a C++20 library for computing Fourier-transformed Green’s functions via iterated Equations of Motion and Lanczos-based resolvent methods. 

The library provides `XPResolvent` for problems with an \(XP\)-structured operator basis and `GeneralResolvent` for more general dynamical and norm matrices. 

`XPResolvent` also supports computing operator amplitudes for eigenoperators in the special \(XP\)-structured case. 

## Requirements

Required dependencies are a C++20 compiler and Eigen. 

Recommended dependencies are OpenMP, `nlohmann/json`, and CMake 3.30 or newer. 

BLAS and LAPACK are optional and may be used through Eigen for faster linear-algebra routines. 

## Basic Usage

Users typically derive their own class from `XPResolvent` or `GeneralResolvent`. 

The derived class must implement `fill_M()`, `fill_matrices()`, and `create_starting_states()`. 

```cpp
struct MyResolvent : public XPResolvent<double> {
  void fill_M() override { /* fill dynamical matrix blocks */ }
  void fill_matrices() override { /* fill dynamical and norm blocks */ }
  void create_starting_states() override { /* define Lanczos starting states */ }

  MyResolvent() : XPResolvent<double>(/* ... */) {}
};
```

Run the Lanczos-based computation with `compute_collective_modes(n)`, where `n` is the number of Lanczos iterations. 

The result is stored as `ResolventDataWrapper` objects containing the computed Lanczos/resolvent data. 

## Main Classes

- `XPResolvent`: optimized implementation for \(XP\)-structured problems with matrix blocks `K_plus`, `K_minus`, and `L`. 
- `GeneralResolvent`: more general implementation using full matrices `M` and `N`. 
- `XPStartingState`: helper container for phase-like and amplitude-like Lanczos starting states. 
- `Resolvent`: lower-level Lanczos implementation used internally by the resolvent classes. 

## Optional Features

Define `MROCK_IEOM_DO_NOT_PARALLELIZE` to disable OpenMP parallelization. 

Define `MROCK_IEOM_PARALLELIZE_BLOCKMATRIX` to parallelize block-matrix operations with OpenMP. 

Define `MROCK_IEOM_NO_NLOHMANN_JSON` if `nlohmann/json` is unavailable. 

## Example

A worked BCS example is available in `tests/ieom_bcs.cpp`. 

The accompanying Python script `tests/bcs_spectral.py` evaluates spectral functions from the computed continued-fraction data. 

## More Information

For full mathematical background, API details, class templates, preprocessor flags, and the complete BCS walkthrough, see `user_guide`. 

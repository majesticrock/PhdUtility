# mrock - umbrella library

mrock is an umbrella C++20 project that bundles three subprojects: iEoM, symbolic_operators, and utility, and provides a common header/interface target and build integration. 

## Requirements

The project requires
- C++20 toolchain
- Eigen
optional/recommended dependencies used by the subprojects are
- OpenMP
- nlohmann/json
- Boost
- CMake 3.30 or newer 
See the subproject READMEs for exact details. 

## Build & install

You can build using the provided Makefile (targets include configure, build, test, install, clean) which wraps CMake, or call CMake directly as described in the top-level CMakeLists. 

An interactive installer script, test_and_install.sh, is provided to help configure component selection, extra include paths, tests, and the install prefix; the script will configure, build, optionally run tests, and optionally install the project. 

## Configuration options

Component toggles (MROCK_BUILD_UTILITY, MROCK_BUILD_SYMBOLIC_OPERATORS, MROCK_BUILD_IEOM) default to ON and can be changed via CMake or the installer script, and you can pass extra include directories or CMAKE_PREFIX_PATH to locate external dependencies. 

The build system uses a default install prefix (../.mrock relative to the source) unless you set CMAKE_INSTALL_PREFIX. 

## Documentation & examples

Detailed user guides for iEoM and symbolic_operators are provided as user_guide.tex.txt and user_guide_symop.tex.txt, and the utility subproject overview is in README_utility.md.txt. 

Example tests and small demos live under each subproject's tests directory and are exercised by the top-level test and install flow. 

## Quick start

- Configure & build via Makefile:
- make configure; make build; make test; make install. 
- Or use the installer script:
- ./test_and_install.sh (follow prompts). 
- Or run CMake directly:
- cmake -S . -B build && cmake --build build --parallel. 

For API details and class-level documentation, consult the subproject user guides referenced above. 
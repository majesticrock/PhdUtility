# mrock - umbrella library

mrock is an umbrella C++20 project that bundles three subprojects: `iEoM`, `symbolic_operators`, and `utility`, and provides a common header/interface target and build integration. 

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

Component toggles (`MROCK_BUILD_UTILITY`, `MROCK_BUILD_SYMBOLIC_OPERATORS`, `MROCK_BUILD_IEOM`) default to `ON` and can be changed via CMake or the installer script, and you can pass extra include directories or CMAKE_PREFIX_PATH to locate external dependencies. 

The build system uses a default install prefix (../.mrock relative to the source) unless you set `CMAKE_INSTALL_PREFIX`. 

## Documentation & examples

User guides for `iEoM` and `symbolic_operators`, as well as the `utility` subproject overview are provided in the respective directory.

Example tests and small demos live under each subproject's tests directory and are exercised by the top-level test and install flow. 

## Quick start

- Configure & build via Makefile:
- make configure; make build; make test; make install. 
- Or use the installer script:
- ./test_and_install.sh (follow prompts). 
- Or run CMake directly:
- cmake -S . -B build && cmake --build build --parallel. 

For API details and class-level documentation, consult the subproject user guides referenced above. 

### Metadata header

A CMake function generates a small metadata header embedding git, build, and system information at configure time.
By default it writes the header (e.g., `generated/include/mrock/info.h`) from a template and exposes variables with the generated paths for CMake and installation.
The header defines a namespaced `info` struct with constexpr string fields such as `GIT_COMMIT_VERSION`, `GIT_COMMIT_NAME`, `GIT_COMMIT_DATE`, `MAKE_DATE`, and `CXX_COMPILER` that you can include and read at runtime.
If git or other probes fail during generation, the function prints a warning and leaves fallback values (e.g., `"unknown"`), allowing the build to continue.
Invoke the generator via `mrock_generate_information_header` at configure time and include the produced header in your code to access `mrock::info::...` identifiers.
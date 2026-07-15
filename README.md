This project contains various general functionality that I require(d) during my PhD.
The structure is the following:

## python:
Contains various python scripts.
The most relevant script are the continued_fraction scripts and the cpp_continued_fraction within the cpp subfolder.
They work in unison to efficiently evaluate the continued fractions occurring in the iEoM formalism.

### Prerequisites
- pybind11 
- pip

### Installation
- Go to `python/cpp`
- Type `pip install .`

To use the other scripts, add the folder to your include path.

## utility
Contains various C++ functionality, ranging from meta-programming to numerics and a config file reader.
Its heart is the iEoM implementation in Numerics/iEoM.

### Prerequisites
- Eigen (for some of the code) https://eigen.tuxfamily.org/index.php?title=Main_Page or https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
- Boost (for some of the code) https://www.boost.org/
- Cmake (for installation) https://cmake.org/

### Installation
From the repository root, run the WSL-compatible build script:

```sh
bash ./scripts/build_and_test.sh
```

This configures, builds, runs the tests, and installs the project into the local `install/` directory under the repository root. The build is portable and does not assume a machine-specific installation prefix or CPU target.

### Usage
This section is header-only. Therefore, you can simply add `~/usr/local` to your include path and use `#include <mrock/utility/Numerics/iEoM/XPResolvent.hpp>`.

## symbolic_operators
Contains the functionality for commuting bosonic and fermionic creation and annihilation operators.
A rather detailed documentation is provided within the folder.

###  Prerequisites
- Boost (iostreams and serialization, if you want to serialize the output)
- Cmake

### Installation
Same as utility.

### Usage
For the header files, this works as utility: Add `~/usr/local` to your include path then include what you need, e.g., `#include <mrock/symbolic_operators/Term.hpp>`.
Afterwards, however, you need to link against the precompiled library! How to do so depends on your compiler.
It is easy to achieve this, if you are using cmake:
- Add `~/usr/local` to your cmake path.
- Use `find_package(mrock REQUIRED)`
- Use `target_include_directories(YourProject PRIVATE ${mrock_INCLUDE_DIRS})`
- Use `target_link_libraries(YourProject PRIVATE ${mrock_LIBRARIES})`
- You are done!

The recipe without cmake is not tested, but should look something like this:
```sh
g++ -I ~/usr/local/include -L ~/usr/local/lib -lmrock_symbolic_operators -o your_binary your_source.cpp
```
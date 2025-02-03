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
- Eigen (for some of the code)
- Boost (for some of the code)
- Cmake (for installation)

### Installation
Go to the root folder of PhdUtility and type `make install`. This builds and installs the entire library (including symbolic_operators) and installs it to `~/usr/local/`.
The header files are copied to `include` and the precompiled sources to `lib`.

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
- Use `find_package(mrock REQUIRED)` or `find_dependency(mrock REQUIRED)`
- Use `target_include_directories(YourProject PRIVATE ${mrock_INCLUDE_DIRS})`
- Use `target_link_libraries(YourProject PRIVATE ${mrock_LIBRARIES})`
- You are done!

The recipe without cmake is not tested, but should look something like this:
```sh
g++ -I ~/usr/local/include -L ~/usr/local/lib -lmrock_symbolic_operators -o your_binary your_source.cpp
```
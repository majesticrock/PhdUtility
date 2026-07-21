# Python Tools for Evaluating iEoM Simulation Data

Small Python package for post-processing simulation data generated with the iEoM library.  
It loads resolvent data, evaluates continued fractions, computes spectral densities, classifies bound states, and provides helpers for plotting and filtering full-diagonalization data.

For details, see the accompanying TeX documentation.

## Requirements

- Python
- numpy
- pandas
- setuptools
- wheel
- pybind11
- C++17-capable compiler
- Boost headers

If Boost is installed in a non-standard location, set:

```bash
export BOOST_INCLUDE_DIR=/path/to/boost
```

## Installation

From the `python` directory, run:

```bash
pip install -e .
```

## Data directory

The data directory can be configured with:

```bash
export MROCK_DATA_DIR=/path/to/data
```

## Basic usage

```python
from mrock.get_data import DataLoader, hubbard_params
import mrock.continued_fraction as cf

loader = DataLoader()

data = loader.load_panda(
    model="hubbard",
    subdir="cube/test",
    file="resolvents.json.gz",
    **hubbard_params(T=0.0, U=-2.5, V=-0.1)
)

resolvents = cf.ContinuedFraction(data)
```

## Main components

- `DataLoader`: load compressed JSON and pickle data.
- `ContinuedFraction`: evaluate continued fractions and spectral densities.
- `FullDiagPurger`: clean and classify full-diagonalization data.
- `cpp_continued_fraction`: C++ backend for numerical routines.

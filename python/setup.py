from pathlib import Path
import os
import platform

from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import numpy


ROOT = Path(__file__).parent.resolve()

extra_compile_args = []
extra_include_dirs = []

if platform.system() == "Windows":
    extra_compile_args = ["/O2"]
    
    boost_include_dir = os.environ.get("BOOST_INCLUDE_DIR")
    extra_include_dirs.append(boost_include_dir)
else:
    extra_compile_args = ["-O3"]

    boost_include_dir = os.environ.get("BOOST_INCLUDE_DIR")
    if boost_include_dir:
        extra_include_dirs.append(boost_include_dir)


ext_modules = [
    Pybind11Extension(
        "mrock.cpp_continued_fraction",
        sources=[
            str(ROOT / "sources" / "cpp_continued_fraction.cpp"),
        ],
        include_dirs=[
            numpy.get_include(),
            *extra_include_dirs,
        ],
        cxx_std=17,
        extra_compile_args=extra_compile_args,
    ),
]


setup(
    name="mrock",
    version="1.0.0",
    description="Python and C++ continued-fraction utilities.",
    packages=find_packages(include=["mrock", "mrock.*"]),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=[
        "numpy",
        "matplotlib",
    ],
    zip_safe=False,
)
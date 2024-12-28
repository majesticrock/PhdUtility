from setuptools import setup, Extension
import pybind11
import numpy
import platform

extra_compile_args = []
if platform.system() == "Windows":
    extra_compile_args = ['/O2']
else:
    extra_compile_args = ['-O3']

ext_modules = [
    Extension(
        'cpp_continued_fraction',  # name of the module
        ['cpp_continued_fraction.cpp'],  # source files
        include_dirs=[pybind11.get_include(), numpy.get_include()],
        language='c++',
        extra_compile_args=extra_compile_args
    ),
]

setup(
    name='cpp_continued_fraction',
    ext_modules=ext_modules,
    zip_safe=False,
)

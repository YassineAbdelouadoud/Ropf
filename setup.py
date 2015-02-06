from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import sys
from distutils.core import setup
from Cython.Build import cythonize

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='Ropf',
    version='0.1',
    packages=find_packages(),
    tests_require=['pytest'],
    cmdclass={'test': PyTest},
    url='',
    license='',
    include_package_data=True,
    author='Yassine Abdelouadoud',
    author_email='yassine.abdelouadoud@gmail.com',
    description='Optimal Power Flow in Radial Networks',
    install_requires=["numpy", "cvxpy", "Cython"],
    ext_modules=cythonize(
                "powerflow.pyx",            # Cython source
                sources=["network.cpp"],  # additional source file(s)
                language="c++",             # generate C++ code
                )
)

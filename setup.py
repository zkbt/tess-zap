#!/usr/bin/env python

# The template for this setup.py came from Tim Morton,
# who I understand took it from Dan F-M. And then Geert
# Barentsen and Christina Hedges helped explain a few
# more neat tips. Thanks all!


import os
import sys
from setuptools import setup, find_packages

# Prepare and send a new release to PyPI
if "release" in sys.argv[-1]:
    os.system("python setup.py sdist")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/tesszap*")
    sys.exit()


# return the README as a string
def readme():
    with open('README.md') as f:
        return f.read()

# a little kludge to be able to get the version number from the package
import sys
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__TESSZAPSETUP__ = True
import tesszap
version = tesszap.__version__

setup(name = "tesszap",
    version = version,
    description = "Simple tool to mitigate cosmic rays from astronomical timeseries.",
    long_description = readme(),
    author = "Zachory K. Berta-Thompson",
    author_email = "zach.bertathompson@colorado.edu",
    url = "https://github.com/zkbt/tess-zap",
    packages = find_packages(),
    package_data = {'tesszap':[]},
    include_package_data=False,
    scripts = [],
    classifiers=[
      'Intended Audience :: Science/Research',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    install_requires=['numpy', 'matplotlib'],
    zip_safe=False,
    license='MIT',
)

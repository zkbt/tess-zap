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
    # uncomment this to test out on test.pypi.com/project/tess-zap
    # os.system("twine upload --repository-url https://test.pypi.org/legacy/ dist/*")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/tesszap*")
    sys.exit()

# a little kludge to be able to get the version number from the package
import sys
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__TESSZAPSETUP__ = True
import tesszap
version = tesszap.__version__

setup(name = "tess-zap",
    version = version,
    description = "Tool to simulate the effects of the TESS cosmic ray mitigation on light curves.",
    long_description = "For usage, installation, and discussion, please visit https://github.com/zkbt/tess-zap",
    author = "Zach Berta-Thompson",
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
    install_requires=['numpy', 'matplotlib', 'scipy', 'astropy', 'tqdm'],
    zip_safe=False,
    license='MIT',
)

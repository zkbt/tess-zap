__version__ = '0.0.1'

# specify whether we're calling this from within setup.py
try:
    __TESSZAPSETUP__
except NameError:
    __TESSZAPSETUP__ = False

if not __TESSZAPSETUP__:
    # (run this stuff if it's not form within setup.py)
    pass

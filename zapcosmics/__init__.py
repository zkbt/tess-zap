__version__ = '0.0'

# specify whether we're calling this from within setup.py
try:
    __ZAPCOSMICS_SETUP__
except NameError:
    __ZAPCOSMICS_SETUP__ = False

if not __ZAPCOSMICS_SETUP__:
    # (run this stuff if it's not form within setup.py)
    pass

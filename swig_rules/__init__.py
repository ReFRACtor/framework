# Just import any files we find in this directory, rather than listing
# everything

from ._swig_wrap import *

# Pull this into all.py. Not sure if we want to do this always or not,
# but this at least gives a way for external inclusion to only include
# a portion of this.
if False:
    import os as _os
    import glob as _glob

    for _i in _glob.glob(_os.path.dirname(__file__) + "/*.py"):
        exec('from .' + _os.path.basename(_i).split('.')[0] + ' import *')

    del _os
    del _glob
    del _i


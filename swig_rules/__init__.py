# Just import any files we find in this directory, rather than listing
# everything

from ._swig_wrap import *

# Pull this into all.py. Not sure if we want to do this always or not,
# but this at least gives a way for external inclusion to only include
# a portion of this.

# Geocal is set up to assume everything is in the top module. Not sure how
# to sort this out, but for now just continue doing this. It isn't clear
# how useful the "all.py" is anyways.
#if False:
if True:
    import os as _os
    import glob as _glob

    for _i in _glob.glob(_os.path.dirname(__file__) + "/*.py"):
        exec('from .' + _os.path.basename(_i).split('.')[0] + ' import *')

    del _os
    del _glob
    del _i


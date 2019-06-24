# Just import any files we find in this directory, rather than listing
# everything

from __future__ import absolute_import
import os as _os
import glob as _glob

from ._swig_wrap import *

for _i in _glob.glob(_os.path.dirname(__file__) + "/*.py"):
    exec('from .' + _os.path.basename(_i).split('.')[0] + ' import *')

del _os
del _glob
del _i


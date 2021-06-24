# For now, we just fail if refractor_framework_swig isn't available. 
# We really pretty much need this, if it isn't installed there isn't
# much of anything we can do. If there was any reason to change this in
# the future we could put some kind of conditional in here.
from refractor_framework_swig import *
have_refractor_framework_swig = True

import os as _os
import re as _re
import glob as _glob

for _i in _glob.glob(_os.path.dirname(__file__) + "/*.py"):
    mname = _os.path.basename(_i).split('.')[0]
    # Don't load ipython, which is ipython magic extensions, or unit tests
    if(not mname == "ipython" and
       not mname == "cython_try" and
       not _re.search('_test', mname)):
        exec("from .%s import *" % mname)

from .config import *
from .director import *
from .executor import *
from .input import *
from .output import *
        
del _i
del _re
del _os
del _glob

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

del _i
del _re
del _os
del _glob

import os as _os
import glob as _glob

for _i in _glob.glob(_os.path.dirname(__file__) + "/*.py"):
    mname = _os.path.basename(_i).split('.')[0]
    if(mname != "all"):
        exec(f"from .{mname} import *")

del _os
del _glob
del _i

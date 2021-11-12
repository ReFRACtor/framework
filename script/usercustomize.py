# This is a usercustomize.py file that will process the ".pth" files that
# might be in the directory (e.g., easy-install.pth created by pip).
# We set these up to be preferred to any system paths, so code in this
# directory gets loaded even if there is a different version available in
# the path.
#
# The issue here is that we can point to our install directory using
# PYTHONPATH, but the ".pth" files don't get processed.
import site
import os
import sys

before_len = len(sys.path)
site.addsitedir(os.path.dirname(__file__))
# Want the new paths to actually appear first, we want to prefer these
# to whatever else might be in the system
if(len(sys.path) != before_len):
    t = sys.path[before_len:]
    t.extend(sys.path[:before_len])
    sys.path = t
if False:    
    print('Paths:')
    for p in sys.path:
        print('  ', p)
    print


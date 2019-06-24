from __future__ import absolute_import
import sys
import glob
import os
import re

# Depending on how we are built, we might or might not have SWIG available.
# The python modules need to work either way, so we have this simple test
# and load if found.
if False:
    try:
        import refractor_swig as framework
        have_refractor_swig = True
    except ImportError:
        have_refractor_swig = False
else:
# The way we use this now, pretty much need swig to do anything. Better to
# just let any error (e.g. link error) go through
    import refractor_swig as framework
    have_refractor_swig = True

# Make sure we can safely import matplotlib without getting an error
# (see this module for details on this)
from . import safe_matplotlib_import

# Normal __all__ doesn't work here, because we also need to include the
# refractor_swig stuff (with __all__, *only* the listed modules are
# included). So we just explicitly import the modules we want here.

# Don't automatically import these modules, they may use C interface
# stuff and should not be available unless directly imported
NO_AUTO_IMPORT = ["__init__",]

for i in glob.glob(os.path.dirname(__file__) + "/*.py"):
    mname = os.path.basename(i).split('.')[0]
    # Don't automatically import test routines
    if(not re.match('.*_test', mname)) and (not mname in NO_AUTO_IMPORT):
        exec('from .%s import *' % mname)

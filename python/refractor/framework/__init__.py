# For now, we just fail if refractor_framework_swig isn't available. 
# We really pretty much need this, if it isn't installed there isn't
# much of anything we can do. If there was any reason to change this in
# the future we could put some kind of conditional in here.
from refractor_framework_swig import *
have_refractor_framework_swig = True

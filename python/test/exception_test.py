from test_support import *

def test_exception():
    t = rf.FpException("test")
    # We added a back trace, so the exception doesn't return a constant
    # value
    #assert t.what() == "test"
        

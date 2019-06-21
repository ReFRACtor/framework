from test_support import *

def test_exception():
    t = rf.FpException("test")
    assert t.what() == "test"
        

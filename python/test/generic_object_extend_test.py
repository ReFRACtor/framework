from refractor.generic_object_extend import *
from test_support import *
from refractor_swig import *

def test_subobj(sample_forward_model):
    #print(list(sample_forward_model.subobject()))
    assert len(list(sample_forward_model.find_all_subobject_of_type(IlsInstrument))) == 1
    print(sample_forward_model.find_subobject_of_type(IlsInstrument))
    
    

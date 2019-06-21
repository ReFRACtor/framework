from test_support import *
import pickle

def test_psigma_pickle():
    '''Test pickling of PressureSigma class.'''
    psigma = rf.PressureSigma([0,0,0],[0.3,0.6,1.0], 10, True)
    t = pickle.dumps(psigma)
    psigma2 = pickle.loads(t)
    assert_array_almost_equal(psigma.a, psigma2.a)
    assert_array_almost_equal(psigma.b, psigma2.b)

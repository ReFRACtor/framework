from .test_support import *
import pickle

def test_psigma_pickle():
    '''Test pickling of PressureSigma class.'''
    psigma = rf.PressureSigma([0,0,0],[0.3,0.6,1.0], 10)
    t = pickle.dumps(psigma)
    psigma2 = pickle.loads(t)
    assert psigma.a == approx(psigma2.a)
    assert psigma.b == approx(psigma2.b)

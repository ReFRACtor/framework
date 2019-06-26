from .test_support import *

def test_gsl_bard():
    '''Use GSL solver to find solution solve bard problem.'''
    nlls = rf.BardNLLSProblem()
    nlls.parameters = [1.0,1.0,1.0]
    solver =rf.NLLSSolverGSL(nlls, 100, 1e-5, 1e-5, 1e-5, False)
    assert nlls.cost >  0.001
    solver.solve()
    assert nlls.cost >  0.001
    # This illustrated a memory leak message from SWIG that we will want
    # to fix.
    print(solver.accepted_points)

def test_fm(config_forward_model):
    '''This illustrates our lack of support of boost::optional<blitz::Range>'''
    print(config_forward_model.stacked_pixel_range(0))

def test_to_vector():
    '''This illustrates a problem with to_vector.'''
    a = rf.ArrayAd_double_1(2,3)
    a[0] = rf.AutoDerivativeDouble(1)
    a[1] = rf.AutoDerivativeDouble(2)
    assert a.to_vector()[0].value == a[0].value

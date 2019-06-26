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
    

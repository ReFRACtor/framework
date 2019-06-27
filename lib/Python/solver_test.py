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
    '''This illustrates our previous lack of support of 
    boost::optional<blitz::Range> - fixed now'''
    print(config_forward_model.stacked_pixel_range(0))

def test_to_vector():
    '''This illustrates a problem with to_vector.'''
    a = rf.ArrayAd_double_1(2,3)
    a[0] = rf.AutoDerivativeDouble(1)
    a[1] = rf.AutoDerivativeDouble(2)
    # This gives an error, because the vector in a.to_vector() gets freeded
    # while we still have the [0] in python - leading to a pointer to free
    # memory. Not a general problem, most of our vectors are native types
    # (e.g., double) or boost::shared_ptr which doesn't have the same issue.
    # A std::vector<AutoDerivative<double> > is a special case that causes
    # problems.
    #
    # Run with "valgrind --log-file=temp.log --malloc-fill=0xff
    # --free-fill=0xff `which python` `which pytest`" to force freeded data
    # to NaN so we can fully illustrate this problem
 
    if(False):
        assert a.to_vector()[0].value == a[0].value

    # This succeeds, because the vector is still around in "t" while the
    # value is being used.
    t = a.to_vector()
    assert t[0].value == a[0].value
    

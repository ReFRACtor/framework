from .test_support import *

def test_gsl_bard():
    '''Use GSL solver to find solution solve bard problem.'''
    nlls = rf.BardNLLSProblem()
    nlls.parameters = [1.0,1.0,1.0]
    solver =rf.NLLSSolverGSL(nlls, 100, 1e-5, 1e-5, 1e-5, False)
    assert nlls.cost >  0.01
    solver.solve()
    assert nlls.cost <  0.01
    # This illustrated a previous memory leak message from SWIG 
    # - this has been fixed now.
    pt = solver.accepted_points
    assert len(pt) > 3

def test_fm(config_forward_model):
    '''This illustrates our previous lack of support of 
    boost::optional<blitz::Range> - fixed now'''
    assert config_forward_model.stacked_pixel_range(0).first() == 0
    assert config_forward_model.stacked_pixel_range(0).last() == 835

def test_to_list():
    '''This ensures there are no problems with ArrayAd::to_list'''
    a = rf.ArrayAd_double_1(2,3)
    a[0] = rf.AutoDerivativeDouble(1)
    a[1] = rf.AutoDerivativeDouble(2)
    # If a vector was used instead of a list then the following would
    # result in an error, because the vector gets freeded
    # while we still have the [0] in python - leading to a pointer to free
    # memory. Not a general problem, most of our vectors are native types
    # (e.g., double) or boost::shared_ptr which doesn't have the same issue.
    # A std::vector<AutoDerivative<double> > is a special case that causes
    # problems. Therefore, its best to create a list on the Python side
    # iteratively instead of letting SWIG convert a vector to a list.
    #
    # Run with "valgrind --log-file=temp.log --malloc-fill=0xff
    # --free-fill=0xff `which python` `which pytest`" to force freeded data
    # to NaN so we can fully illustrate this problem.
    # This passes now, but failed when a vector was used previously.
 
    assert a.to_list()[0].value == a[0].value

    # Even with a vector, this would succeeds, because the vector is still
    # around in "t" while the value is being used.
    t = a.to_list()
    assert t[0].value == a[0].value
    

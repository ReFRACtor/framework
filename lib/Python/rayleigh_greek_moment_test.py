from .test_support import *

def test_rayleigh_greek_moment():
    '''Test RayleighGreekMoment'''
    expected = np.array([[1,     0, 0, 0, 0, 0],
                         [1e-11, 0, 0, 1.3968144385817844, 0, 0],
                         [0.47936288771635682, 2.8761773262981407, 0, 0, 
                          1.1741944765321404, 0]])
    assert_array_almost_equal(rf.RayleighGreekMoment.array(), expected)
    

from test_support import *
import numpy as np
from scipy import interpolate
import pickle
from pytest import approx

def test_absco_coeff():
    '''Compare Matt's original python test with our code.
    Matt had code at /home/mthill/muses/py-retrieve/mthill/new_absco_tables/sample_interpolation_code_2.py that calculated a few points. We just compare
    our results with his calculation.
    '''
        
    fname = f"{absco_data_dir}/coeff/O3_19840-37879_v0.0_init_new.nc"
    f = rf.AbscoCoeff(fname)
    # Calculate using our C++ code
    def absco_rf_calc(wn, t, p):
        return f.absorption_cross_section(wn, rf.DoubleWithUnit(p, "Pa"), rf.DoubleWithUnit(t, "K"), rf.ArrayWithUnit_double_1([0], "dimensionless")).value

    assert absco_rf_calc(28840.544030722263,220,138.1155) == approx(8.69761749e-23)
    assert absco_rf_calc(28840.544030722263,230,138.1155) == approx(1.10401266e-22)
    assert absco_rf_calc(28840.544030722263,225,138.1155) == approx(9.86887206e-23)
    

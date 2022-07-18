from test_support import *
import numpy as np
from scipy import interpolate
import pickle

@skip
def test_create_matt_data():
    '''We unfortunately can't use both h5py and our framework hdf reader
    in the same program - there is some kind of conflict. Read the
    data from matt's absco and save as a pickle, which we can then use in
    the comparison test.'''
    print("hi there")
    import h5py

    fname = f"{absco_data_dir}/coeff/O3_19840-37879_v0.0_init_new.nc"
    f = h5py.File(fname, 'r')
    # Subset we read of the temperature and pressure. This is pretty much
    # arbitrary, just done to restrict size that we read and also matches
    # Matt's initial test.
    tsub = slice(4,9)
    psub = slice(45,50)
    nrepvecs = f['Cross_Section_repvecs'].shape[0]
    # We may need to modify the interpolation done here, but our results
    # shouldn't be a strong function of how we do the interpolation
    G = interpolate.interp1d(f['Spectral_Grid'][:],
                                   f['Cross_Section_repvecs'][:, :],
                                   kind='linear', axis=1)
    F = [interpolate.interp2d(f['T_layer'][psub, tsub].transpose(),
                              f['P_layer'][psub, tsub].transpose(),
                              f['Cross_Section_coeffs'][i,tsub,psub])
         for i in range(nrepvecs)]
    pickle.dump([G,F],open(f"{unit_test_data}/matt_absco.p","wb"))
    
def test_absco_coeff():
    '''Compare Matt's original python test with our code.

    This depends on data we save in the previous test. I had originally
    thought of just adding this to the repository, but the data is ~1 GB
    which is a bit large. We'll likely only need this test occasionally, 
    so we'll at least for now check to see if the support file is there
    and only run if found.
    '''
    if(not os.path.exists(f"{unit_test_data}/matt_absco.p")):
        raise SkipTest("Need to have matt_absco.p data available to run")
    fname = f"{absco_data_dir}/coeff/O3_19840-37879_v0.0_init_new.nc"
    print(f"{unit_test_data}/matt_absco.p")
    print(fname)
    f = rf.AbscoCoeff(fname)

    G,F = pickle.load(open(f"{unit_test_data}/matt_absco.p","rb"))

    def absco_matt_calc(wn, t, p):
        abscorepvec_interp = G(wn)
        abscocoeff_interp_list = np.array([F[i](t, p) for i in range(len(F))])
        absco_interp = np.matmul(abscorepvec_interp,abscocoeff_interp_list)
        return absco_interp[0]

    # Calculate using our C++ code
    def absco_rf_calc(wn, t, p):
        return f.absorption_cross_section(wn, rf.DoubleWithUnit(p, "Pa"), rf.DoubleWithUnit(t, "K"), rf.ArrayWithUnit_double_1([0], "dimensionless")).value
    
    print(absco_matt_calc(28800.24008156947, 225, 148))
    print(absco_rf_calc(28800.24008156947, 225, 148))
    

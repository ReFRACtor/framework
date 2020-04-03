# Tests the Python interface to the PCA Eigen Solver because 
# this component necessitated added SWIG rules for passing
# a vector of Blitz arrays

from .test_support import *
import os

import pytest

@pytest.fixture(scope='session')
def opt_props():
    opt_prop_filename = os.path.join(unit_test_data, "in/pca/pca_optical_properties.h5")

    opt_props = rf.HdfFile(opt_prop_filename)

    return opt_props

def test_eigen_solver_generic(opt_props):

    data_list = []
    data_list.append(opt_props.read_double_2d("gas_optical_depth"))
    data_list.append(opt_props.read_double_2d("single_scattering_albedo"))

    pca_solver = rf.PCAEigenSolverGeneric(data_list, 4)

def test_bad_element(opt_props):

    data_list = []
    data_list.append(opt_props.read_double_2d("gas_optical_depth"))
    data_list.append("")

    with pytest.raises(TypeError) as exc_info:
        pca_solver = rf.PCAEigenSolverGeneric(data_list, 4)

    assert exc_info.value.args[0] == "iter_to_vector_of_arrays: object is not a numpy object at index: 1"

def test_bad_iterator(opt_props):

    with pytest.raises(TypeError) as exc_info:
        pca_solver = rf.PCAEigenSolverGeneric(None, 4)

    assert exc_info.value.args[0] == "iter_to_vector_of_arrays: passed object is not iterable"

    with pytest.raises(TypeError) as exc_info:
        pca_solver = rf.PCAEigenSolverGeneric(100, 4)

    assert exc_info.value.args[0] == "iter_to_vector_of_arrays: passed object is not iterable"

    with pytest.raises(TypeError) as exc_info:
        pca_solver = rf.PCAEigenSolverGeneric("", 4)

    assert exc_info.value.args[0] == "iter_to_vector_of_arrays: unicode objects are not supported"

def test_incorrect_dimensions(opt_props):

    data_list = []
    data_list.append(np.zeros(0))

    with pytest.raises(TypeError) as exc_info:
        pca_solver = rf.PCAEigenSolverGeneric(data_list, 4)

    assert exc_info.value.args[0] == "iter_to_vector_of_arrays: incorrect dimension of numpy object at index: 0, expected dimension: 2"

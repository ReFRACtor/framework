import os.path
import shutil

from test_support import *

import netCDF4
import numpy as np
from pytest import approx


@pytest.mark.skipif(absco_data_dir is None or
                    not os.path.exists(f"{absco_data_dir}/aer/CH4_06200-06201_v0.0_init.nc"),
                    reason="ABSCO aer test file not found")
def test_values_outside_valid_range_are_read_as_zero(tmp_path):
    '''LBLRTM writes tiny negative optical depths in some line wings, so a table can hold
    negative absorption coefficients. AbscoAer reads values outside the valid_range of the
    table as zero.'''
    src = f"{absco_data_dir}/aer/CH4_06200-06201_v0.0_init.nc"
    fname = str(tmp_path / "CH4_negative.nc")
    shutil.copy(src, fname)
    with netCDF4.Dataset(fname, "r+") as d:
        d.set_auto_mask(False)
        k = d.variables["Cross_Section"]
        grid = d.variables["Spectral_Grid"][:]
        p = d.variables["P_layer"][:]
        t = d.variables["T_layer"][:]
        # a layer and temperature with data, and a grid point near the middle of the table
        iw = len(grid) // 2
        it, il = np.argwhere(np.isfinite(p.T))[0]
        wn, pressure, temperature = float(grid[iw]), float(p[il, it]), float(t[il, it])
        positive = float(k[iw, it, il])
        assert positive > 0
        lo, hi = k.valid_range
        assert lo == 0 and hi > positive
        k[iw, it, il] = -positive
        k[iw + 1, it, il] = 10 * hi

    f = rf.AbscoAer(fname)
    nobroad = rf.ArrayWithUnit_double_1([0], "dimensionless")
    got = f.absorption_cross_section(wn, rf.DoubleWithUnit(pressure * 100.0, "Pa"),
                                     rf.DoubleWithUnit(temperature, "K"), nobroad).value
    assert got == 0.0
    wn_hi = float(grid[iw + 1])
    assert f.absorption_cross_section(wn_hi, rf.DoubleWithUnit(pressure * 100.0, "Pa"),
                                      rf.DoubleWithUnit(temperature, "K"), nobroad).value == 0.0

    # the untouched table still returns the original value
    g = rf.AbscoAer(src)
    assert g.absorption_cross_section(wn, rf.DoubleWithUnit(pressure * 100.0, "Pa"),
                                      rf.DoubleWithUnit(temperature, "K"), nobroad).value == approx(positive)

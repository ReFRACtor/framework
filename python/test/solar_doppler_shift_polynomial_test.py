from test_support import *
import os
import datetime;
import time;

def test_solar_doppler_shift_polynomial():
    # Time as a datetime. Note you can also use a float as a unix time
    # stamp if that is more convenient. We include fractions of a second
    tm = rf.Time.parse_time("2006-09-14T12:27:22.000100Z")

    p = rf.SolarDopplerShiftPolynomial(tm,
                                rf.DoubleWithUnit(77.1828918457,rf.Unit("deg")),
                                rf.DoubleWithUnit(74.128288269,rf.Unit("deg")),
                                rf.DoubleWithUnit(167.495071411,rf.Unit("deg")),
                                rf.DoubleWithUnit(416,rf.Unit("m")))

    assert p.solar_distance.value == approx(1.0060305651331354)

    wn = rf.SpectralDomain([12929.94, 12979.93, 13029.93, 13079.93, 
                         13129.93, 13179.93])
    expected = np.array([12929.919173650407,
                         12979.909093131146,
                         13029.909012595777,
                         13079.908932060409,
                         13129.908851525041,
                         13179.908770989672])
    assert p.doppler_stretch(wn).data == approx(expected)
    

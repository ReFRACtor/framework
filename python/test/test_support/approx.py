from refractor.framework import *
import pytest
from _pytest.python_api import ApproxBase

def _double_with_unit_rep(self):
    return str(self)

DoubleWithUnit.__repr__ = _double_with_unit_rep

class ApproxDoubleWithUnit(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.value ==
                pytest.approx(actual.value, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok)
                and self.expected.units == actual.units)

class ApproxSpectralDomain(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.data ==
                pytest.approx(actual.data, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.sample_index ==
                pytest.approx(actual.sample_index, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.units == actual.units)

class ApproxSpectralRange(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.data ==
                pytest.approx(actual.data, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.uncertainty ==
                pytest.approx(actual.uncertainty, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.units == actual.units)
    
class ApproxArrayWithUnit(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.value ==
                pytest.approx(actual.value, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.units == actual.units)

class ApproxArrayAd(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.value ==
                pytest.approx(actual.value, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.jacobian ==
                pytest.approx(actual.jacobian, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok))

class ApproxArrayAdWithUnit(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.value ==
                approx(actual.value, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.units == actual.units)

class ApproxAutoDerivative(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.value ==
                pytest.approx(actual.value, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.gradient ==
                pytest.approx(actual.gradient, rel=self.rel,
                              abs = self.abs, nan_ok = self.nan_ok))

class ApproxAutoDerivativeWithUnit(ApproxBase):
    def __repr__(self):
        return "approx(%s)" % self.expected
    
    def __eq__(self, actual):
        return (self.expected.value ==
                approx(actual.value, rel=self.rel,
                       abs = self.abs, nan_ok = self.nan_ok) and
                self.expected.units == actual.units)
    
    
# We copy the handling that pytest does to add a "approx" function
def approx(expected, rel=None, abs=None, nan_ok=False):
    if(isinstance(expected, DoubleWithUnit)):
        cls = ApproxDoubleWithUnit
    elif(isinstance(expected, SpectralDomain)):
        cls = ApproxSpectralDomain
    elif(isinstance(expected, SpectralRange)):
        cls = ApproxSpectralRange
    elif(isinstance(expected, ArrayWithUnit_double_1) or
         isinstance(expected, ArrayWithUnit_double_2) or
         isinstance(expected, ArrayWithUnit_double_3) or
         isinstance(expected, ArrayWithUnit_double_4)):
        cls = ApproxArrayWithUnit
    elif(isinstance(expected, ArrayAd_double_1) or
         isinstance(expected, ArrayAd_double_2) or
         isinstance(expected, ArrayAd_double_3) or
         isinstance(expected, ArrayAd_double_4)):
        cls = ApproxArrayAd
    elif(isinstance(expected, ArrayAdWithUnit_double_1) or
         isinstance(expected, ArrayAdWithUnit_double_2) or
         isinstance(expected, ArrayAdWithUnit_double_3) or
         isinstance(expected, ArrayAdWithUnit_double_4)):
        cls = ApproxArrayAdWithUnit
    elif(isinstance(expected, AutoDerivativeDouble)):
        cls = ApproxAutoDerivative
    elif(isinstance(expected, AutoDerivativeWithUnitDouble)):
        cls = ApproxAutoDerivativeWithUnit
    else:
        # Try pytest.approx
        return pytest.approx(expected, rel=rel, abs=abs, nan_ok=nan_ok)

    return cls(expected, rel=rel, abs=abs, nan_ok=nan_ok)

__all__ = ["approx"]                  

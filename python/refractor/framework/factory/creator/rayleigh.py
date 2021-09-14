import os

import numpy as np

from .base import Creator
from .. import param

import refractor.framework as rf

class RayleighYoung(Creator):
    
    pressure = param.InstanceOf(rf.Pressure)
    altitude = param.ObjectVector("altitude")
    constants = param.InstanceOf(rf.Constant)

    def create(self, **kwargs):

        pressure = self.pressure()
        altitude = self.altitude()
        constants = self.constants()
        
        return rf.RayleighYoung(pressure, altitude, constants)

class RayleighBodhaine(Creator):
    
    pressure = param.InstanceOf(rf.Pressure)
    altitude = param.ObjectVector("altitude")
    constants = param.InstanceOf(rf.Constant)

    def create(self, **kwargs):

        pressure = self.pressure()
        altitude = self.altitude()
        constants = self.constants()
        
        return rf.RayleighBodhaine(pressure, altitude, constants)

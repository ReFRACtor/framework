import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class StateVector(Creator):

    def create(self, **kwargs):
        return rf.StateVector()

import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class StateVector(Creator):

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

        self.register_to_receive(rf.SubStateVectorObserver)

    def receive(self, rt_obj):
        pass

    def create(self, **kwargs):
        return rf.StateVector()

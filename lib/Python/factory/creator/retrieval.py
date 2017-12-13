import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class StateVector(Creator):

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

        self.register_to_receive(rf.SubStateVectorObserver)

    def receive(self, rt_obj):
        print(rt_obj.__class__.__name__, hasattr(rt_obj, "sub_state_identifier"))

    def create(self, **kwargs):
        return rf.StateVector()

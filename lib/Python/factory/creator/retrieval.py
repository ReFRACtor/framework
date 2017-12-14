import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class StateVector(Creator):

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

        self.register_to_receive(rf.SubStateVectorObserver)

        self.sv_observers = []

    def receive(self, rt_obj):

        if not hasattr(rt_obj, "sub_state_identifier"):
            raise AttributeError("The recieved SubStateVector class: %s does not have a sub_state_identifier attribute" % rt_obj.__class__.__name__)

        if rt_obj.coefficient.value.shape[0] > 0 and rt_obj not in self.sv_observers:
            self.sv_observers.append(rt_obj)

    def create(self, **kwargs):
        sv = rf.StateVector()

        # Register observers into state vector and gather 
        ig_values = []
        for observer in self.sv_observers:
            sv.add_observer(observer)
            ig_values.append(observer.coefficient.value[np.nonzero(sv_obs.used_flag_value)])

        ig = np.concatenate(ig_values)

        if ig.shape[0] != sv.observer_claimed_size:
            raise ValueError("The initial guess vector size %d does not match expected state vector size %d" % (ig.shape[0], sv.observer_claimed_size))
        
        sv.update_state(ig)

        return sv



from refractor.executor import StrategyExecutor
from refractor.output.base import OutputBase
from refractor_swig import ObserverPtrNamedSpectrum

class CaptureRadiance(ObserverPtrNamedSpectrum, OutputBase):

    def __init__(self):
        # Required to initialize director
        ObserverPtrNamedSpectrum.__init__(self)

        self.convolved_spectrum = []
        self.high_res_spectrum = []

    def notify_update(self, named_spectrum):
        if named_spectrum.name == "convolved":
            self.convolved_spectrum.append(named_spectrum.copy())
        elif named_spectrum.name == "high_res_rt":
            self.high_res_spectrum.append(named_spectrum.copy())

class ComparisonExecutor(StrategyExecutor):
    "Executor that will store produced radiances internally for later comparison against offline values"

    def __init__(self, config_filename, **kwargs):
        super().__init__(config_filename, **kwargs)

        self.captured_radiances = CaptureRadiance()

    def attach_output(self, config_inst, step_index=0, simulation=False):
        super().attach_output(config_inst, step_index, simulation)

        config_inst.forward_model.add_observer_and_keep_reference(self.captured_radiances)

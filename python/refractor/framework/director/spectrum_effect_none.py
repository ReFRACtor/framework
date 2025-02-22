import refractor.framework as rf
from refractor.framework import SpectrumEffectImpBase

class SpectrumEffectNone(SpectrumEffectImpBase):
    """
    Class to hold a SpectrumEffect which passes through original Spectrum

    SpectrumEffectList is a list of vectors of SpectrumEffects.
    Each of the vectors in the list is size num_channels.
    For SpectrumEffects that are not applicable to all channels, this effect can be used for
    non-applicable channels.
    """

    def apply_effect(self, spec, forward_model_grid):
        pass

    def name(self):
        return "passthrough"

    def clone(self):
        return type(self)()

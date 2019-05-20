import numpy as np

from .base import Creator
from .value import CreatorFlaggedValue
from .. import param

from refractor import framework as rf

class PerChannelMixin(object):

    num_channels = param.Scalar(int)

    def check_num_channels(self, check_num):
        "Issue a warning when the number of channels in some data source does not match the expected number"

        if not self.num_channels() == check_num:
            raise param.ParamError("Number of channels for data source: %d does not match expected number: %d" % (check_num, self.num_channels()))

class SpectralWindowRange(Creator):
    
    window_ranges = param.ArrayWithUnit(dims=3)
    bad_sample_mask = param.Array(dims=2, required=False)
    
    def create(self, **kwargs):

        win_ranges = self.window_ranges()
        bad_sample_mask = self.bad_sample_mask()

        if bad_sample_mask is not None:
            return rf.SpectralWindowRange(win_ranges, bad_sample_mask)
        else:
            return rf.SpectralWindowRange(win_ranges)

class MicroWindowRanges(Creator):
    """Creates a SpectralWindowRange object from a list of microwindows. The microwindows can 
    can be in any order, but they need to be contained in the full_ranges."""

    # List of full spectral ranges for each band
    # n_channels x 2
    full_ranges = param.ArrayWithUnit(dims=2)

    # List of microwindows to be mapped to spectral channels
    # n_microwindows x 2
    micro_windows = param.ArrayWithUnit(dims=2, required=False)

    def create(self, **kwargs):

        full_ranges = self.full_ranges()
        micro_windows = self.micro_windows()
        num_channels = full_ranges.rows

        # If micro_windows are not supplied use the full_ranges as the
        # window values.
        if micro_windows is None:
            micro_windows = full_ranges

        # Convert microwindows to the units
        mw_conv = micro_windows.convert_wave(full_ranges.units)

        # Map microwindows to spectral channels
        mapping = { chan_idx: [] for chan_idx in range(num_channels) }
        max_num_mw = 0
        for mw_range in mw_conv.value:
            # convert_wave may have reversed order
            mw_beg, mw_end = sorted(mw_range)

            match_found = False
            for chan_idx, (full_beg, full_end) in enumerate(full_ranges.value):
                if mw_beg >= full_beg and mw_end <= full_end:
                    mapping[chan_idx].append( (mw_beg, mw_end) )
                    match_found = True
                    max_num_mw = max(max_num_mw, len(mapping[chan_idx]))
                    break

            if not match_found:
                raise Exception('Could not find channel bounding microwindow: ({}, {})'.format(mw_beg, mw_end))

        # Create microwindow ranges
        win_ranges = np.zeros((num_channels, max_num_mw, 2))
        for chan_idx in range(num_channels):
            for mw_idx, mw_range in enumerate(mapping[chan_idx]):
                win_ranges[chan_idx, mw_idx, :] = mw_range

        # Assign to configuration
        return rf.SpectralWindowRange(rf.ArrayWithUnit_double_3(win_ranges, full_ranges.units))

class SpectrumSamplingBase(Creator, PerChannelMixin):

    high_res_spacing = param.Choice(param.ArrayWithUnit(dims=1), param.DoubleWithUnit())

    def spacing(self):
        "Returns the array with unit value for spacing specified regardless if it was specified as an array or scalar value with unit"

        spacing_val = self.high_res_spacing()
        num_channels = self.num_channels()

        if isinstance(spacing_val, rf.DoubleWithUnit):
            # Create an ArrayWithDouble matching the number of channels used
            spacing_used = rf.ArrayWithUnit_double_1(np.full(num_channels, spacing_val.value), spacing_val.units)
        else:
            self.check_num_channels(spacing_val.value.shape[0])
            spacing_used = spacing_val

        return spacing_used

class FixedSpacingSpectrumSampling(SpectrumSamplingBase):

    def create(self, **kwargs):
        return rf.SpectrumSamplingFixedSpacing(self.spacing())

class UniformSpectrumSampling(SpectrumSamplingBase):

    def create(self, **kwargs):
        return rf.UniformSpectrumSampling(self.spacing())

class NonuniformSpectrumSampling(SpectrumSamplingBase):

    channel_domains = param.Iterable()

    def create(self, **kwargs):
        domains = self.channel_domains()
        full_spec_spacing = rf.SpectrumSamplingFixedSpacing(self.spacing())

        if len(domains) != self.num_channels():
            raise param.ParamError("Number of channel domains %d does not match the number of channels %d" % (len(domains), self.num_channels()))

        for idx, dom in enumerate(domains):
            if not isinstance(dom, rf.SpectralDomain):
                raise param.ParamError("Channel domain value at index %d is not a instance of SpectralDomain" % idx)

        return rf.NonuniformSpectrumSampling(*domains, full_spec_spacing)

class SpectrumEffectList(Creator):

    effects = param.Iterable(str)
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        all_effects = []
        for effect_name in self.effects():
            self.register_parameter(effect_name, param.Iterable())
            all_effects.append(self.param(effect_name))

        # Map these into an outer vector for each channel, with an inner vector for each effect
        spec_eff = rf.vector_vector_spectrum_effect()
        for chan_index in range(self.num_channels()):
            per_channel_eff = rf.vector_spectrum_effect()

            for effect in all_effects:
                per_channel_eff.push_back(effect[chan_index])

            spec_eff.push_back(per_channel_eff)

        return spec_eff

class FluorescenceEffect(CreatorFlaggedValue):
    
      reference_point = param.ArrayWithUnit(dims=1)
      cov_unit = param.InstanceOf(rf.Unit)
      which_channels = param.Iterable()

      atmosphere = param.InstanceOf(rf.AtmosphereOco)
      observation_zenith = param.ArrayWithUnit(dims=1)
      stokes_coefficient = param.Array(dims=2)
      num_channels = param.Scalar(int)

      def create(self, **kwargs):

          atm = self.atmosphere()
          lza = self.observation_zenith()
          ref_point = self.reference_point()[0]
          coeff = self.value()
          used_flag = self.retrieval_flag()
          stokes_coeff = rf.StokesCoefficientConstant(self.stokes_coefficient())
          cov_unit = self.cov_unit()
          num_channels = self.num_channels()

          # Create array for all channels, but filled with null objects (None)
          fluoresence = [None] * num_channels

          for chan_idx in self.which_channels():
              chan_idx = int(chan_idx)

              if chan_idx >= num_channels:
                  raise ParamError("Channel index {} exceeds number of channels {}".format(chan_idx, num_channels))

              fluoresence[chan_idx] = rf.FluorescenceEffect(coeff, used_flag, atm, stokes_coeff, lza[chan_idx], chan_idx, ref_point, cov_unit)

          return fluoresence

class ForwardModel(Creator):

    instrument = param.InstanceOf(rf.Instrument)
    spec_win = param.InstanceOf(rf.SpectralWindow)
    radiative_transfer = param.InstanceOf(rf.RadiativeTransfer)
    spectrum_sampling = param.InstanceOf(rf.SpectrumSampling)
    spectrum_effect = param.ObjectVector("vector_spectrum_effect")

    def create(self, **kwargs):
        fm = rf.StandardForwardModel(self.instrument(),
                                     self.spec_win(), 
                                     self.radiative_transfer(),
                                     self.spectrum_sampling(), 
                                     self.spectrum_effect())

        fm.setup_grid()

        return fm

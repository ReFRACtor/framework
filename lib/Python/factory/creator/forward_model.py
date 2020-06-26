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
    # Defaults to full_ranges if not supplied or equal to None
    micro_windows = param.Choice(param.ArrayWithUnit(dims=2), param.NoneValue(), required=False)

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
        previous_effects = {}
        for effect_name in self.effects():
            self.register_parameter(effect_name, param.Iterable())
            per_chan_effects = self.param(effect_name, **previous_effects)
            all_effects.append(per_chan_effects)
            # Store effects as they are encountered so that subsequent steps can utilize previous ones
            previous_effects[effect_name] = per_chan_effects

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

      atmosphere = param.InstanceOf(rf.AtmosphereStandard)
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

class RamanSiorisEffect(CreatorFlaggedValue):

    albedo = param.Array(dims=1)

    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    relative_azimuth = param.ArrayWithUnit(dims=1)
    atmosphere = param.InstanceOf(rf.AtmosphereStandard)
    solar_model = param.Choice(param.ObjectVector("vector_spectrum_effect"), param.Iterable(rf.SolarModel))
    num_channels = param.Scalar(int)

    do_upwelling = param.Scalar(bool, default=True)
    padding_fraction = param.Scalar(float, default=0.10)
    jacobian_perturbation = param.Scalar(float, default=0.001)

    def create(self, solar_model=None, **kwargs):

        scale_factor = self.value()
        used_flag = self.retrieval_flag()

        albedo = self.albedo()

        solar_zenith = self.solar_zenith()
        obs_zenith = self.observation_zenith()
        rel_azimuth = self.relative_azimuth()
        atmosphere = self.atmosphere()

        # Get value from keyword if not passed through
        if solar_model is None:
            solar_model = self.solar_model()

        padding_fraction = self.padding_fraction()
        do_upwelling = self.do_upwelling()
        jac_perturb = self.jacobian_perturbation()

        raman_effect = []
        for chan_index in range(self.num_channels()):
            raman_effect.append( rf.RamanSiorisEffect(scale_factor[chan_index], bool(used_flag[chan_index]), chan_index,
                                 solar_zenith[chan_index], obs_zenith[chan_index], rel_azimuth[chan_index],
                                 atmosphere, solar_model[chan_index], albedo[chan_index], 
                                 padding_fraction, do_upwelling, jac_perturb) )

        return raman_effect

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

class OssForwardModel(Creator):
    absorber = param.Iterable(rf.AbsorberVmr)
    pressure = param.InstanceOf(rf.Pressure)
    temperature = param.InstanceOf(rf.Temperature)
    surface_temperature = param.InstanceOf(rf.SurfaceTemperature)
    ground = param.InstanceOf(rf.GroundPiecewise)
    # TODO: Discuss below with respect to channels
    # observation_zenith = param.DoubleWithUnit()
    # solar_zenith = param.DoubleWithUnit()
    # latitude = param.DoubleWithUnit()
    # altitude = param.DoubleWithUnit()
    observation_zenith = param.ArrayWithUnit()
    solar_zenith = param.ArrayWithUnit()
    latitude = param.ArrayWithUnit()
    altitude = param.ArrayWithUnit()
    # TODO: some configs set this up at vector of altitude instead of ArrayWithUnit from L1B
    # altitude = param.ObjectVector("altitude")
    lambertian = param.Scalar(bool, default=True)
    sel_file = param.Scalar(str)
    od_file = param.Scalar(str)
    sol_file = param.Scalar(str)
    fix_file = param.Scalar(str)
    ch_sel_file = param.Scalar(str)
    max_chans =  param.Scalar(int, default=20000)

    ## TEMP ##
    # Temporary to attach a ForwardModelSpectralGrid until changes are made
    # to fix the requirement of a ForwardModelSpectralGrid object to create a retrieval
    # TODO: Remove this
    instrument = param.InstanceOf(rf.Instrument)
    spec_win = param.InstanceOf(rf.SpectralWindow)
    spectrum_sampling = param.InstanceOf(rf.SpectrumSampling)
    ## TEMP ##

    def create(self, **kwargs):

        vmrs = rf.vector_absorber_vmr()
        for vmr_obj in self.absorber():
            vmrs.push_back(vmr_obj)

        # TODO: discuss support for multiple channels/bands
        chan_idx = 0
        obs_zen = rf.DoubleWithUnit(self.observation_zenith().value[chan_idx], self.observation_zenith().units)
        sol_zen = rf.DoubleWithUnit(self.solar_zenith().value[chan_idx], self.solar_zenith().units)
        latitude = rf.DoubleWithUnit(self.latitude().value[chan_idx], self.latitude().units)
        surf_alt = rf.DoubleWithUnit(self.altitude().value[chan_idx], self.altitude().units)
        fm = rf.OssForwardModel(vmrs, self.pressure(), self.temperature(),
                                self.surface_temperature(), self.ground(),
                                obs_zen, sol_zen, latitude, surf_alt, self.lambertian(),
                                self.sel_file(), self.od_file(), self.sol_file(), self.fix_file(),
                                self.ch_sel_file(), self.max_chans())

        fm.setup_grid()

        ## TEMP ##
        # Temporary to attach a ForwardModelSpectralGrid until changes are made
        # to fix the requirement of a ForwardModelSpectralGrid object to create a retrieval
        # TODO: Remove this
        fm.spectral_grid = rf.ForwardModelSpectralGrid(self.instrument(), self.spec_win(), self.spectrum_sampling())
        ## TEMP ##

        return fm

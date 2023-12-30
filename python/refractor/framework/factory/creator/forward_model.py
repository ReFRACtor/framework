from bisect import bisect

import numpy as np

from .base import Creator
from .types import RetrievalComponents
from .. import param

import refractor.framework as rf

class PerChannelMixin(object):

    num_channels = param.Scalar(int)

    def check_num_channels(self, check_num):
        "Issue a warning when the number of channels in some data source does not match the expected number"

        if not self.num_channels() == check_num:
            raise param.ParamError("Number of channels for data source: %d does not match expected number: %d" % (check_num, self.num_channels()))

class SpectralWindowRange(Creator):
    
    window_ranges = param.ArrayWithUnit(dims=3)
    # Mask of instrument samples that should not be included in spectral regions
    bad_sample_mask = param.Choice(param.Array(dims=2),
                                   param.Iterable(np.ndarray), required=False)
    
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

    # Mask of instrument samples that should not be included in spectral regions
    bad_sample_mask = param.Choice(param.Array(dims=2),
                                   param.Iterable(np.ndarray), required=False)

    def create(self, **kwargs):

        full_ranges = self.full_ranges()
        micro_windows = self.micro_windows()
        bad_sample_mask = self.bad_sample_mask()

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
        win_ranges = rf.ArrayWithUnit(win_ranges, full_ranges.units)

        if bad_sample_mask is not None:
            return rf.SpectralWindowRange(win_ranges, bad_sample_mask)
        else:
            return rf.SpectralWindowRange(win_ranges)

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
        spec_eff = []
        for chan_index in range(self.num_channels()):
            per_channel_eff = []

            for effect in all_effects:
                per_channel_eff.append(effect[chan_index])

            spec_eff.append(per_channel_eff)

        return spec_eff

class FluorescenceEffect(Creator):
    
      coefficients = param.Array(dims=1)
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
          coeff = self.coefficients()
          stokes_coeff = rf.StokesCoefficientConstant(self.stokes_coefficient())
          cov_unit = self.cov_unit()
          num_channels = self.num_channels()

          # Create array for all channels, but filled with null objects (None)
          fluoresence = [None] * num_channels

          for chan_idx in self.which_channels():
              chan_idx = int(chan_idx)

              if chan_idx >= num_channels:
                  raise ParamError("Channel index {} exceeds number of channels {}".format(chan_idx, num_channels))

              fluoresence[chan_idx] = rf.FluorescenceEffect(coeff, atm, stokes_coeff, lza[chan_idx], chan_idx, ref_point, cov_unit)

          return fluoresence

class RamanSiorisEffect(Creator):
    instrument = param.InstanceOf(rf.Instrument)
    spec_win = param.InstanceOf(rf.SpectralWindow)
    spectrum_sampling = param.InstanceOf(rf.SpectrumSampling)

    scale_factors = param.Array(dims=1)
    solar_zenith = param.ArrayWithUnit(dims=1)
    observation_zenith = param.ArrayWithUnit(dims=1)
    relative_azimuth = param.ArrayWithUnit(dims=1)
    atmosphere = param.InstanceOf(rf.AtmosphereStandard)
    solar_model = param.Choice(param.ObjectVector("vector_spectrum_effect"), param.Iterable(rf.SolarModel))
    num_channels = param.Scalar(int)

    do_upwelling = param.Scalar(bool, default=True)

    def create(self, solar_model=None, **kwargs):
        # TODO: Better way to pick up high resolution grid?
        spectral_grid = rf.ForwardModelSpectralGrid(self.instrument(), self.spec_win(), self.spectrum_sampling())

        scale_factor = self.scale_factors()

        solar_zenith = self.solar_zenith()
        obs_zenith = self.observation_zenith()
        rel_azimuth = self.relative_azimuth()
        atmosphere = self.atmosphere()

        # Get value from keyword if not passed through
        if solar_model is None:
            solar_model = self.solar_model()

        do_upwelling = self.do_upwelling()

        mapping = rf.StateMappingLinear()

        raman_effect = []
        for chan_index in range(self.num_channels()):
            raman_effect.append( rf.RamanSiorisEffect(spectral_grid.high_resolution_grid(chan_index),
                                 scale_factor[chan_index], chan_index,
                                 solar_zenith[chan_index], obs_zenith[chan_index], rel_azimuth[chan_index],
                                 atmosphere, solar_model[chan_index], mapping, do_upwelling) )

        return raman_effect

class Cloud3dEffect(Creator):
    offset = param.Array(dims=1)
    slope = param.Array(dims=1)

    desc_band_name = param.Iterable()
    num_channels = param.Scalar(int)

    def create(self, **kwargs):

        offset = self.offset()
        slope = self.slope()
        band_name = self.desc_band_name()

        cloud_3d_effect = []
        for chan_index in range(self.num_channels()):
            cloud_3d_effect.append( rf.Cloud3dEffect(offset[chan_index], slope[chan_index], band_name[chan_index]) )

        return cloud_3d_effect

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
    #max_chans =  param.Scalar(int, default=20000)

    ## TEMP ##
    # Temporary to attach a ForwardModelSpectralGrid until changes are made
    # to fix the requirement of a ForwardModelSpectralGrid object to create a retrieval
    # TODO: Remove this
    instrument = param.InstanceOf(rf.Instrument)
    spec_win = param.InstanceOf(rf.SpectralWindow)
    spectrum_sampling = param.InstanceOf(rf.SpectrumSampling)
    ## TEMP ##

    def __init__(self, *vargs, **kwargs):
        super().__init__(*vargs, **kwargs)

        self.register_to_receive(RetrievalComponents, deregister_after_create=False)

        # Expected order of retrieval components so that they match the order that OssForwardModel provides
        self.ret_object_order = [ rf.Temperature,
                                  rf.SurfaceTemperature,
                                  rf.AbsorberVmr,
                                  rf.GroundPiecewise ]

        self.fm = None
        self.receive_is_registered = False

    def receive(self, rec_obj):
        if not self.receive_is_registered:
            self.setup_retrieval_levels(list(rec_obj.values()))
            self.receive_is_registered = True

    def check_retrieval_order(self, retrieval_components):

        # Go through and assign an index to each of the retrieval components based on the above ordering
        configured_order = []
        for ret_obj in retrieval_components:
            found_type = False
            for idx, obj_type in enumerate(self.ret_object_order):
                if isinstance(ret_obj, obj_type):
                    configured_order.append(idx)
                    found_type = True
                    break

            if not found_type:
                raise param.ParamError(f"Retrieval component configured that OssForwardModel can not handle: {ret_obj.__class__.__name__}")

        # Check that the configured order is in increasing order
        last_idx = -1
        for pos_idx, config_idx in enumerate(configured_order):
            if config_idx < last_idx:
                ret_obj = retrieval_components[pos_idx]
                raise param.ParamError(f"Retrieval component {ret_obj.__class__.__name__} is in the wrong order.")
            last_idx = config_idx

    def temperature_retrieval_levels(self, retrieval_components, fm_press):

        temperature = self.temperature()

        if temperature in retrieval_components:
            temp_press = temperature.pressure_profile()

            temp_ret_levels = []
            for ret_idx, press in enumerate(temp_press):
                # TODO: New way to lookup
                # if temperature.used_flag_value[ret_idx]:
                    temp_ret_levels.append(bisect(fm_press, press))

            return np.array(temp_ret_levels, dtype=int)
        else:
            return np.zeros(0, dtype=int)

    def surface_temperature_flags(self, retrieval_components):

        surface_temperature = self.surface_temperature()

        if surface_temperature in retrieval_components:
            # TODO: New way to lookup
            # return np.array(surface_temperature.used_flag_value, dtype=bool)
            return np.ones(self.fm.num_channels, dtype=bool)
        else:
            # Jacobian for all sensor channels disabled
            return np.zeros(self.fm.num_channels, dtype=bool)

    def gas_retrieval_levels(self, retrieval_components, fm_press):

        absorber = self.absorber()

        # Gases
        gas_ret_levels = []
        for vmr_obj in self.absorber():
            vmr_press = vmr_obj.coeff_pressure_profile
            # TODO: New way to lookup
            """
            if vmr_obj in retrieval_components and np.any(vmr_obj.used_flag_value):

                vmr_levels = []
                for ret_idx, press in enumerate(vmr_press):
                    if vmr_obj.used_flag_value[ret_idx]:
                        vmr_levels.append( bisect(fm_press, press) )

                gas_ret_levels.append( np.array(vmr_levels, dtype=int) )
            else:
                gas_ret_levels.append( np.zeros(0, dtype=int) )
            """
            if vmr_obj in retrieval_components:
                vmr_levels = []
                for ret_idx, press in enumerate(vmr_press):
                    vmr_levels.append( bisect(fm_press, press) )
                gas_ret_levels.append( np.array(vmr_levels, dtype=int) )
            else:
                gas_ret_levels.append( np.zeros(0, dtype=int) )
        return gas_ret_levels

    def ground_flags(self, retrieval_components):

        ground = self.ground()

        if ground in retrieval_components:
            if isinstance(ground, rf.GroundEmissivityPiecewise):
                # TODO: New way to lookup
                # emiss_flags = np.where(ground.used_flag_value)[0]
                emiss_flags = np.ones(ground.spectral_points().rows, dtype=int)
            else:
                emiss_flags = np.zeros(0, dtype=int)

            if isinstance(ground, rf.GroundLambertianPiecewise):
                # TODO: New way to lookup
                # refl_flags = np.where(ground.used_flag_value)[0]
                refl_flags = np.ones(ground.spectral_points().rows, dtype=int)
            else:
                refl_flags = np.zeros(0, dtype=int)

        return emiss_flags, refl_flags

    def setup_retrieval_levels(self, retrieval_components):

        pressure = self.pressure()

        fm_press = pressure.pressure_grid().value.value

        self.check_retrieval_order(retrieval_components)

        temp_ret_levels = self.temperature_retrieval_levels(retrieval_components, fm_press)
        surf_temp_flags = self.surface_temperature_flags(retrieval_components)
        gas_ret_levels = self.gas_retrieval_levels(retrieval_components, fm_press)
        emiss_flags, refl_flags = self.ground_flags(retrieval_components)

        # Create retrieval flags object and set up forward model object
        ret_flags = rf.OssRetrievalFlags(temp_ret_levels, surf_temp_flags, gas_ret_levels, emiss_flags, refl_flags)

        self.fm.setup_retrieval(ret_flags)
        self.fm.setup_grid()

    def create(self, **kwargs):

        vmrs = []
        for vmr_obj in self.absorber():
            vmrs.append(vmr_obj)

        # TODO: discuss support for multiple channels/bands
        chan_idx = 0
        obs_zen = rf.DoubleWithUnit(self.observation_zenith().value[chan_idx], self.observation_zenith().units)
        sol_zen = rf.DoubleWithUnit(self.solar_zenith().value[chan_idx], self.solar_zenith().units)
        latitude = rf.DoubleWithUnit(self.latitude().value[chan_idx], self.latitude().units)
        surf_alt = rf.DoubleWithUnit(self.altitude().value[chan_idx], self.altitude().units)
        self.fm = rf.OssForwardModel(vmrs, self.pressure(), self.temperature(),
                                     self.surface_temperature(), self.ground(),
                                     obs_zen, sol_zen, latitude, surf_alt, self.lambertian(),
                                     self.sel_file(), self.od_file(), self.sol_file(), self.fix_file(),
                                     self.ch_sel_file())


        ## TEMP ##
        # Temporary to attach a ForwardModelSpectralGrid until changes are made
        # to fix the requirement of a ForwardModelSpectralGrid object to create a retrieval
        # TODO: Remove this
        self.fm.spectral_grid = rf.ForwardModelSpectralGrid(self.instrument(), self.spec_win(), self.spectrum_sampling())
        ## TEMP ##

        return self.fm

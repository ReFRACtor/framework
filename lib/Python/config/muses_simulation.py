import os
import re
from copy import deepcopy
from glob import glob

import numpy as np

from netCDF4 import Dataset
from scipy.io import readsav

import refractor.factory.creator as creator
from refractor import framework as rf

from refractor.muses_osp import MUSES_File

class MusesSimConfig(object):

    # Must be defined in inheriting instances
    grid_units = None

    def __init__(self, config_def, num_channels):

        self.config_def = config_def
        self.num_channels = num_channels

        # Make sure configuration routine is only called once
        self.configured = False

    def microwindows_from_grid(self, sim_grid):
        
        # Create microwindows from simulation_grid
        dom_diff = (sim_grid[1:] - sim_grid[:-1])
        wend = np.where(dom_diff > dom_diff[0])[0]
        wbeg = np.concatenate(([0], wend+1))
        wend = np.concatenate((wend, [sim_grid.shape[0]-1]))

        mw_points = []
        for mw_b, mw_e in zip(wbeg, wend):
            mw_points.append( (sim_grid[mw_b], sim_grid[mw_e]) )

        return mw_points

    def set_microwindows(self, micro_windows):

        # Use the correct creator
        self.config_def['spec_win'].update({
            'creator': creator.forward_model.MicroWindowRanges,
            'micro_windows': rf.ArrayWithUnit(np.array(micro_windows), self.grid_units),
        })

    def disable_retrieval(self):

        # Disable retrieval
        self.config_def['order'].remove('retrieval')

        # Disable L1B input
        self.config_def['input'] = {}

    def setup_scenario(self, latitude=0, longitude=0, obs_azimuth=0, obs_zenith=0, solar_azimuth=0, solar_zenith=1e-6, surface_height=0, sample_grid=None):

        self.config_def['scenario'] =  {
            'creator': creator.base.SaveToCommon,
            "latitude": rf.ArrayWithUnit(np.full(self.num_channels, latitude, dtype=float), "deg"),
            "longitude": rf.ArrayWithUnit(np.full(self.num_channels, longitude, dtype=float), "deg"),
            "observation_azimuth": rf.ArrayWithUnit(np.full(self.num_channels, obs_azimuth, dtype=float), "deg"),
            "observation_zenith": rf.ArrayWithUnit(np.full(self.num_channels, obs_zenith, dtype=float), "deg"),
            "solar_azimuth": rf.ArrayWithUnit(np.full(self.num_channels, solar_azimuth, dtype=float), "deg"),
            "solar_zenith": rf.ArrayWithUnit(np.full(self.num_channels, solar_zenith, dtype=float), "deg"),
            "stokes_coefficient": np.tile(np.array([1.0, 0.0, 0.0, 0.0]), (self.num_channels, 1)),
            "surface_height": rf.ArrayWithUnit(np.full(self.num_channels, surface_height, dtype=float), "m"),
            "sample_grid": "NOT_YET_SET",
        }

    def setup_sample_grid(self, sample_grid):

        self.config_def['scenario']['sample_grid'] = sample_grid

        # Use sample grid from scenario
        self.config_def['instrument']['dispersion'] = {
            'creator': creator.instrument.SampleGridSpectralDomain,
            'value': {
                'creator': creator.value.NamedCommonValue,
                'name': 'sample_grid',
            },
        }

    def setup_atmosphere_profiles(self, atm_gas_list, profile_getter):

        # Set pressure grid, reverse and convert mb -> Pa
        self.config_def['atmosphere']['pressure']['pressure_levels'] =  profile_getter("pressure") * 100
        self.config_def['atmosphere']['pressure']['value'] = np.array([profile_getter("pressure")[-1] * 100])

        # Set temperature values
        self.config_def['atmosphere']['temperature']['temperature_levels'] = profile_getter("TATM")

        # Set absorbers

        # Get names of gases we actually know about
        absco_gas_list = [ fn.split("_")[0] for fn in os.listdir(os.environ['ABSCO_PATH']) ]

        # Ignore due to errors that occur when including these gases, possible errors in ABSCO or the profile
        ignore_gas_list = ['CCL4', 'PAN']

        used_gas_list = []
        for in_gas_name in atm_gas_list:
            absco_gas_name = in_gas_name

            absco_gas_name = re.sub('^HDO$', 'H2O', absco_gas_name)

            if not absco_gas_name in absco_gas_list:
                print("No ABSCO available for gas named: {}".format(absco_gas_name))
                continue
            elif absco_gas_name in ignore_gas_list:
                print("Ignoring troublesome gas: {}".format(absco_gas_name))
                continue

            if not absco_gas_name in used_gas_list:
                used_gas_list.append(absco_gas_name)

            if absco_gas_name in self.config_def['atmosphere']['absorber']:
                gas_def = self.config_def['atmosphere']['absorber'][absco_gas_name]
                if gas_def['vmr']['value'] is not None:
                    print("Adding {} to existing gas {}".format(in_gas_name, absco_gas_name))
                    gas_def['vmr']['value'] += profile_getter(in_gas_name)
                else:
                    gas_def['vmr']['value'] = profile_getter(in_gas_name)
            else:
                # Must make a deepcopy so we are not giving each gas a reference to the default definition that changing
                # would make all have the same values
                gas_def = deepcopy(self.config_def['atmosphere']['absorber']['default_gas_definition'])
                gas_def['vmr']['value'] = profile_getter(in_gas_name)
                self.config_def['atmosphere']['absorber'][absco_gas_name] = gas_def

        self.config_def['atmosphere']['absorber']['gases'] = used_gas_list

    def setup_emissivity(self, emiss_grid, emiss_values):

        # Setup emissivity values using the picewise method
        self.config_def['atmosphere']['ground']['grid'] = rf.ArrayWithUnit(emiss_grid, "cm^-1")
        self.config_def['atmosphere']['ground']['value'] = emiss_values

    def setup_surface_temperature(self, surface_temp):

        self.config_def['atmosphere']['surface_temperature']['value'] = rf.ArrayWithUnit(np.full(self.num_channels, surface_temp, dtype=float), "K")

    def setup_surface_albedo(self, albedo_values):

        self.config_def['atmosphere']['ground']['lambertian']['value'] = albedo_values

    def configure_scenario(self, sample_grid):
        # Instruments must define the way they configure the scenario, they can use setup_scenario
        # to set values
        raise NotImplementedError("configure_scenario must be implemented in an inheriting class")

    def configure_micro_windows(self):
        raise NotImplementedError("configure_micro_windows must be implemented in an inheriting class")

    def configure_convolved_grid(self):
        raise NotImplementedError("configure_convolved must be implemented in an inheriting class")

    def configure_atmosphere(self):
        raise NotImplementedError("configure_atmosphere must be implemented in an inheriting class")

    def configure(self):

        if self.configured:
            raise Exception("Configuration setup already called")

        # Disable retrieval
        self.disable_retrieval()

        self.configure_scenario()

        self.configure_micro_windows()

        self.configure_sample_grid()

        self.configure_atmosphere()

        self.configured = True

class MusesCaseDirectoryConfig(MusesSimConfig):
    
    def __init__(self, muses_case_dir, config_def, num_channels, atm_gas_list=None):
        super().__init__(config_def, num_channels)

        self.muses_case_dir = muses_case_dir

        self.atm_gas_list = atm_gas_list

    def configure_micro_windows(self):

        radiance_fn = os.path.join(self.muses_case_dir, "Products/Products_Radiance-FM-pantest.nc")
        with Dataset(radiance_fn, "r") as rad_file_data:
            sim_grid = rad_file_data['/FREQUENCY'][:]

        # Set microwindows from simulation grid
        microwindows = self.microwindows_from_grid(sim_grid)

        self.set_microwindows(microwindows)

    def configure_sample_grid(self):

        radiance_fn = os.path.join(self.muses_case_dir, "Products/Products_Radiance-FM-pantest.nc")
        with Dataset(radiance_fn, "r") as rad_file_data:
            sim_grid = rad_file_data['/FREQUENCY'][:]

        sample_grid = [ rf.SpectralDomain(sim_grid, rf.Unit(self.grid_units)) ]
        while len(sample_grid) < self.num_channels:
            # Pad out extra unused channels
            sample_grid.append(rf.SpectralDomain(np.array([]), rf.Unit(self.grid_units)))

        self.setup_sample_grid(sample_grid)

    def configure_atmosphere(self):

        # Load atmosphere values
        atm_fn = glob(os.path.join(self.muses_case_dir, "Input/Initial_State/State_AtmProfiles_*.asc"))[0]
        atm_profiles = MUSES_File(atm_fn, as_struct=True)

        column_names_lower = [ c.lower() for c in atm_profiles ]

        def atmosphere_profile(col_name):
            col_idx = column_names_lower.index(col_name)
            return atm_profiles[atm_profiles.column_names[col_idx]]

        if self.atm_gas_list is None:
            # First 3 elements are Level, Pressure, TATM, after that gas names
            atm_gas_list = atm_profiles.column_names[3:]
        else:
            atm_gas_list = self.atm_gas_list

        self.setup_atmosphere_profiles(atm_gas_list, atmosphere_column)

class MusesThermalCaseDirectoryConfig(MusesSimConfig):

    def configure(self):

        super().configure()

        # Load surface data
        surface_file = os.path.join(self.muses_case_dir, "Input/Initial_State/Emissivity1.asc")
        surface_data = MUSES_File(surface_file, as_struct=True)

        emiss_grid = surface_data.data['Frequency']
        emiss_values = surface_data.data['Emissivity']

        self.setup_emissivity(emiss_grid, emiss_values)

        self.setup_surface_temperature(surface_data.header['Surface_Temperature'])

class MusesUipSimConfig(MusesSimConfig):

    def __init__(self, muses_uip_file, config_def, num_channels, atm_gas_list=None):
        super().__init__(config_def, num_channels)

        self.muses_uip_file = muses_uip_file

        self.muses_inputs = readsav(muses_uip_file, python_dict=True)
        self.uip = self.muses_inputs['uip']

        self.atm_gas_list = atm_gas_list

    def configure_micro_windows(self):

        # Set up microwindow ranges
        mw_all = self.uip['microwindows_all'][0]

        microwindows = []
        for start, end, spacing in zip(mw_all['start'], mw_all['endd'], mw_all['monospacing']):
            microwindows.append( (start, end) )

        self.set_microwindows(microwindows)

    def configure_sample_grid(self):

        # Set up microwindow ranges
        mw_all = self.uip['microwindows_all'][0]

        sample_grid = []
        for start, end, spacing in zip(mw_all['start'], mw_all['endd'], mw_all['monospacing']):
            sample_grid.append( rf.SpectralDomain(np.arange(start, end, spacing), rf.Unit("cm^-1")) )

        # Set up scenario
        while len(sample_grid) < self.num_channels:
            # Pad out extra unused channels
            sample_grid.append(rf.SpectralDomain(np.array([]), rf.Unit("cm^-1")))

        self.setup_sample_grid(sample_grid)

    def configure_atmosphere(self):

        # Helper for pulling information out of the atmosphere matrix
        def atmosphere_column(param_name):
            param_list = [ n.decode('UTF-8').lower() for n in self.uip['atmosphere_params'][0] ]
            param_index = param_list.index(param_name.lower())
            return self.uip['atmosphere'][0][:, param_index][::-1]

        if self.atm_gas_list is None:
            # First 2 elements are Pressure and TATM, after that gas names
            gas_list = [ n.decode('UTF-8') for n in self.uip['species'][0] ]
        else:
            gas_list = self.atm_gas_list

        self.setup_atmosphere_profiles(gas_list, atmosphere_column)

class MusesThermalUipSimConfig(MusesUipSimConfig):

    def configure(self):
        super().configure()

        # Load surface data
        emiss_grid = self.uip['emissivity'][0]['frequency'][0]
        emiss_values = self.uip['emissivity'][0]['value'][0]

        self.setup_emissivity(emiss_grid, emiss_values)

        # Surface temperature
        self.setup_surface_temperature(self.uip['surface_temperature'][0])

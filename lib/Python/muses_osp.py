import os
import re
import struct
from functools import lru_cache
from glob import glob
from bisect import bisect_right
import datetime as dt
import warnings

import numpy as np

class MUSES_File(object):

    def __init__(self, filename, as_struct=False, header_only=False):

        self.filename = filename
        self.as_struct = as_struct
        self.header_only = header_only

        self.header = {}
        self.column_names = []
        self.unit_names = []
        self.data = None

        with open(self.filename, "rb") as mf:
            self._read_file(mf)

    def _read_file(self, muses_file):

        self._read_header(muses_file)

        if self.header_only:
            return

        if 'Binary_Endian' in self.header:
            self._read_binary_data(muses_file)
        else:
            self._read_column_names(muses_file)
            self._read_ascii_data(muses_file)

    def _read_header(self, muses_file):

        # Parse heade contents until finding end of header
        curr_line = b""
        while not re.search(rb"^\s*End_of_Header", curr_line):
            curr_line = muses_file.readline()
            
            header_value = re.search(rb"^([\w-]+)\s*=\s*(.*)$", curr_line)
            if header_value:
                key, value = [ s.decode('UTF-8') for s in header_value.groups() ]

                # Scrub any quote marks and extra spaces
                value = value.strip('"').strip()

                # Parse datasize into an array of numbers
                data_size_re = r"\s+x\s*"
                if re.search(data_size_re, value):
                    value = [ int(v) for v in re.split(data_size_re, value) ]

                self.header[key] = value

    def _read_column_names(self, muses_file):

        # Read column names and unit names from lines after end of header
        self.column_names = [ v.decode('UTF-8') for v in muses_file.readline().strip().split() ]
        self.unit_names = [ v.decode('UTF-8') for v in muses_file.readline().strip().split() ]
 
    def _read_ascii_data(self, muses_file, skip_lines=0):

        # Skip column unit labels, then read data
        data = np.genfromtxt(muses_file, names=self.column_names, skip_header=skip_lines)
        
        if self.as_struct:
            self.data = data
        else:
            # Convert from structured array
            if len(data.shape) == 0:
                struct_rows = 1
            else:
                struct_rows = data.shape[0]

            mat = np.empty((struct_rows, len(data.dtype.names)), data.dtype[0])
            for col_idx, name in enumerate(data.dtype.names):
                mat[:, col_idx] = data[name]

            self.data = mat

    def _read_binary_data(self, muses_file):

        data_size = self.header['Data_Size']

        # Assuming data is doubles
        num_doubles = np.prod(data_size)
        num_bytes = 8 * num_doubles

        raw_bytes = muses_file.read(num_bytes)

        if self.header["Binary_Endian"] == "big":
            # big endian
            format_code = ">{}d".format(num_doubles)
        else:
            # little endian
            format_code = "<{}d".format(num_doubles)

        raw_doubles = struct.unpack(format_code, raw_bytes)
        self.data = np.array(raw_doubles).reshape(data_size)

class PressureOSP(object):
    def __init__(self, base_dir, pressure_cutoff=100000, **kwargs):

        self.base_dir = base_dir
        self.pressure_cutoff = pressure_cutoff
 
        # This is really a fixed value for a set of OSP files, but could be changed should a new
        # set of OSPs be created
        self.num_fm_levels = 66

    @property
    def fm_pressure_filename(self):
        "Filename used for accessing the forward model pressure levels"

        press_file = os.path.join(self.base_dir, f"Strategy_Tables/Defaults/TES_baseline_{self.num_fm_levels}.asc")

        if not os.path.exists(press_file):
            raise Exception(f"Could not find forward model pressure file: {press_file}")

        return press_file

    @lru_cache()
    def _fm_pressure_all(self):
        "Forward model pressure grid"

        press_file = self.fm_pressure_filename

        # Reverse to pressure increasing order, convert hPa -> Pa
        mf = MUSES_File(press_file)
        press_data = np.flip(mf.data[:, 0] * 100)

        return press_data

    @property
    def fm_pressure_grid(self):
        "Forward model pressure grid"

        press_data = self._fm_pressure_all()

        if self.pressure_cutoff is not None:
            press_data = press_data[np.where(press_data <= self.pressure_cutoff)]

        return press_data

class SpeciesOSP(PressureOSP):

    def __init__(self, species, base_dir, latitude, longitude, obs_time, instrument_name=None, cov_dir="Covariance", **kwargs):
        super().__init__(base_dir, **kwargs)

        self.species = species
        self.latitude = latitude
        self.longitude = longitude

        if not isinstance(obs_time, dt.datetime):
            raise Exception("obs_time must be a datetime instance")

        self.obs_time = obs_time

        self.cov_dir = cov_dir

        self.instrument_name = instrument_name

        self.nh3_version = "CLN"

    @property
    def climatology_base_dir(self):
        clim_base_dir = os.path.join(self.base_dir, "Climatology", self.species)
        if not os.path.exists(clim_base_dir):
            raise Exception("No climatology species directory found: {} for species".format(clim_base_dir, self.species))
        return clim_base_dir

    @property
    def strategy_base_dir(self):
        if self.instrument_name is None:
            raise Exception("No instrument name supplied to constructor, can not use Strategy_Tables based values.")

        strat_base_dir = os.path.join(self.base_dir, "Strategy_Tables", f"OSP-{self.instrument_name}")
        if not os.path.exists(strat_base_dir):
            raise Exception(f"No Strategy Tables directory found at: {strat_base_dir} for instrument name: {self.instrument_name}")
        return strat_base_dir
 
    @property
    def species_pressure_filename(self):
        "Filename used for accessing the pressure grid associated with a species"

        # Find the associated pressure file
        press_file = os.path.join(self.climatology_base_dir, "Clim_Spec_{}".format(self.species))

        if not os.path.exists(press_file):
            raise Exception("Could not find pressure file: {} for species: {}".format(press_file, self.species))

        return press_file

    @lru_cache()
    def _species_pressure_all(self):

        press_file = self.species_pressure_filename
    
        # Reverse to pressure increasing order, convert hPa -> Pa
        mf = MUSES_File(press_file)
        press_data = np.flip(mf.data[:, 0] * 100)

        return press_data

    def _species_grid_indexes(self):

        press_data = self._species_pressure_all()

        if self.pressure_cutoff is not None:
            grid_indexes = np.where(press_data <= self.pressure_cutoff)
        else:
            grid_indexes = np.arange(press_data.shape[0])

        return grid_indexes

    @property
    def species_pressure_grid(self):
        "The pressure grid associated with a species converted into the format acceptable for ReFRACtor"

        press_data = self._species_pressure_all()

        if self.pressure_cutoff is not None:
            press_data = press_data[np.where(press_data <= self.pressure_cutoff)]

        return press_data[self._species_grid_indexes()]

    def _pick_long_dir(self, longitude_dirs):
        if len(longitude_dirs) == 1:
            return longitude_dirs[0]

        # Convert longitude directory names into a mapping of the directory name to range of latitude values
        long_ranges = { d: [ float(v.strip("E")) for v in d.split("_")] for d in longitude_dirs }

        use_long_dir = None
        for long_dir, long_rang in long_ranges.items():
            long_360 = self.longitude % 360
            if long_360 >= long_rang[0] and long_360 <= long_rang[1]:
                use_long_dir = long_dir
                break

        if use_long_dir is None:
            raise Exception("Could not determine longitude directoryecies: {}".format(self.species))

        return use_long_dir

    def _pick_lat_file(self, latitude_files):

        if len(latitude_files) == 1:
            return latitude_files[0]

        # Convert latitudes with N (North) and S (South) markings to positive and negative latitudes
        def from_ns(lat_str):
            num_val = re.sub(r'[NS]$', '', lat_str)
            if re.search(r'S', lat_str):
                return -float(num_val)
            else:
                return float(num_val)

        # Find the file that exists within the latitude range
        use_lat_file = None
        for lat_file in latitude_files:
            if re.search(r"NH3", lat_file):
                matches = re.search(r'.*(\d\d[NS])_(\d\d[NS])_{}$'.format(self.nh3_version), lat_file)
            else:
                matches = re.search(r'.*(\d\d[NS])_(\d\d[NS])$', lat_file)

            if matches is None:
                continue

            beg, end = [ from_ns(lat_str) for lat_str in matches.groups() ]
            
            if self.latitude >= beg and self.latitude <= end:
                use_lat_file = lat_file
                break

        return use_lat_file

    @property
    def _month_directory_name(self):
        return self.obs_time.strftime("%b").upper()
   
    def _time_directory_name(self):

        species_base_dir = self.climatology_base_dir

        all_year_dirs = sorted([ int(y) for y in filter(lambda d: re.match(r'\d{4}', d), os.listdir(species_base_dir)) ])

        if len(all_year_dirs) > 0:
            year_idx = max(0, bisect_right(all_year_dirs, self.obs_time.year)-1)
            year_dir = str(all_year_dirs[year_idx])
            species_base_dir = os.path.join(species_base_dir, year_dir)
        else:
            year_dir = ''

        all_month_dirs = list(filter(lambda d: re.match(r'[A-Z]{3}', d), os.listdir(species_base_dir)))

        if len(all_month_dirs) == 1:
            month_dir = all_month_dirs[0]
        else:
            best_month_dir = self._month_directory_name
            if best_month_dir not in all_month_dirs:
                raise Exception("Month directory for current data {} not available at {} for species".format(best_month_dir, species_base_dir, self.species))

            month_dir = best_month_dir

        time_dir = "0000UT_2400UT"

        return os.path.join(year_dir, month_dir, time_dir)

    def _temporal_clim_filename(self):
        
        species_base_dir = self.climatology_base_dir

        time_rel_dir = self._time_directory_name()

        time_base_dir = os.path.join(species_base_dir, time_rel_dir)

        if not os.path.exists(time_base_dir):
            raise Exception("No climatology temporal directory found: {} for species: {}".format(time_base_dir, self.species))

        long_dir = self._pick_long_dir(os.listdir(time_base_dir))
        lat_file = self._pick_lat_file(os.listdir(os.path.join(time_base_dir, long_dir)))

        clim_file = os.path.join(time_base_dir, long_dir, lat_file)

        return clim_file

    def _nh3_clim_filename(self):
        return glob(os.path.join(self.climatology_base_dir, self.nh3_version, "Profile_NH3_{}.asc".format(self.nh3_version)))[0]

    @property
    def climatology_filename(self):
        "Filename with the species climatology values matching the object's time and location"

        if self.species == "NH3":
            return self._nh3_clim_filename()
        else:
            return self._temporal_clim_filename()

    @property
    def species_climatology(self):
        "Species climatology data matching the object's time and location converted to the format handled by ReFRACtor"

        clim_file = self.climatology_filename

        mf = MUSES_File(clim_file)

        clim_data = mf.data[:, 0]

        # Convert to pressure increasing order
        clim_data = np.flip(clim_data)

        # Only trim by pressure those OSPs that are of the same size
        # as their associated pressure grid
        if clim_data.shape[0] == self._species_pressure_all().shape[0]:
            return clim_data[self._species_grid_indexes()]
        else:
            return clim_data

    @property
    def fm_climatology(self):

        fm_press = self.fm_pressure_grid
        species_press = self.species_pressure_grid
        species_value = self.species_climatology

        return np.exp(np.interp(np.log(fm_press), np.log(species_press), np.log(species_value)))

    @property
    def retrieval_strategy_filename(self):
        ret_strat_fn = os.path.join(self.strategy_base_dir, f"Species-{self.num_fm_levels}", f"{self.species}.asc")
        if not os.path.exists(ret_strat_fn):
            return None
        return ret_strat_fn

    @property
    def retrieval_strategy_info(self):

        if not hasattr(self, "_ret_strat"):
            if self.retrieval_strategy_filename is None:
                self._ret_strat = {}
            else:
                mf = MUSES_File(self.retrieval_strategy_filename, header_only=True)
                self._ret_strat = mf.header

        return self._ret_strat

    @property
    def use_log(self):

        if self.retrieval_strategy_info.get("mapType", "") == "Log":
            return True
        else:
            return False

    @property
    def retrieval_levels(self):

        ret_lev = self.retrieval_strategy_info.get("retrievalLevels", 0)

        if ret_lev == 0:
            return None
        else:
            fm_press_all = self._fm_pressure_all()
            num_fm_levels = fm_press_all.shape[0]

            # Convert to an array and make zero based indexing
            ret_lev = np.array([ int(l)-1 for l in ret_lev.split(",") ])

            # Reverse order and indexes to match pressure increasing order of ReFRACtor
            ret_lev = np.flip(num_fm_levels - ret_lev - 1)

            if self.pressure_cutoff is not None:
                # Incorporate the pressure cutoff 
                ret_flags = np.zeros(num_fm_levels, dtype=bool)
                ret_flags[ret_lev] = True

                cut_ret_lev = np.where(ret_flags[np.where(fm_press_all <= self.pressure_cutoff)])[0]
            else:
                cut_ret_lev = ret_lev

            return cut_ret_lev

    def _pick_constraint_levels_file(self, filenames):

        ret_levs = self.retrieval_levels

        if ret_levs is not None:
            num_ret_levels = len(ret_levs)
        else:
            num_ret_levels = None

        if num_ret_levels is not None:
            # If defined find the file with the defined number of constraint levels
            for fn in filenames:
                if re.search(r"_{}.asc".format(num_ret_levels), fn):
                    return fn

            raise Exception(f"Could not find a constraint file with {num_ret_levels} retrieval levels")
        else:
            return None

    def _co_star_filename(self, co_dir_suffix, constraint=False):
        "Finds filenames for Covariance and Constraint files which have a similar pattern"

        cov_base_dir = os.path.join(self.base_dir, co_dir_suffix)

        if not os.path.exists(cov_base_dir):
            raise Exception("{} base directory does not exist: {} for species: {}".format(co_dir_suffix, cbase_dir, self.species))

        if constraint:
            glob_pattern = "Constraint_Matrix_{species}_*"
        else:
            glob_pattern = "Covariance_Matrix_{species}_*"

        species_files = glob(os.path.join(cov_base_dir, glob_pattern.format(species=self.species)))

        if len(species_files) == 0:
            return None

        if self.use_log:
            filt_files = filter(lambda f: re.search(r"_Log_", f, re.IGNORECASE), species_files)
        else:
            filt_files = filter(lambda f: re.search(r"_Linear_", f, re.IGNORECASE), species_files)

        if constraint:
            co_file = self._pick_constraint_levels_file(list(filt_files))
        else:
            co_file = self._pick_lat_file(list(filt_files))

        return co_file
 
    @property
    def covariance_filename(self):

        return self._co_star_filename(os.path.join("Covariance", self.cov_dir))

    @property
    def constraint_filename(self):
        
        constraint_fn = self.retrieval_strategy_info.get("constraintFilename", None)

        if constraint_fn is not None:
            # Replace base path with path to our filees
            constraint_fn = constraint_fn.replace("../OSP", self.base_dir)

            num_ret_levels = len(self.retrieval_levels)

            # Replace _87 in the filename with the number of retrieval levels like the MUSES code does
            if self.species == "NH3":
                constraint_fn = constraint_fn.replace("_87", f"_{num_ret_levels}_{self.nh3_version}")
            else:
                constraint_fn = constraint_fn.replace("_87", f"_{num_ret_levels}")

            return constraint_fn
        else:
            return self._co_star_filename(os.path.join("Constraint"), constraint=True)

    @lru_cache()
    def _co_star_read(self, co_file):

        if co_file is None:
            raise Exception("No {} constraint/covariance file found for species: {}".format(self.use_log and "log" or "linear", self.species))

        mf = MUSES_File(co_file)

        # Find pressures column, convert hPa -> Pa, put in increasing pressure order
        press_column = mf.column_names.index("Pressure")
        press_data = np.flip(mf.data[:, press_column] * 100)

        # Extract out actual map
        species_col_names = filter(lambda cn: re.search(self.species, cn), mf.column_names)
        species_col_idx = [ mf.column_names.index(nm) for nm in species_col_names ]
        co_data = mf.data[:, np.array(species_col_idx)]

        # Convert to pressure increasing order
        co_data = np.flip(co_data)

        return press_data, co_data

    @property
    def covariance_pressure(self):

        cov_filename = self.covariance_filename
        cov_press, cov_data = self._co_star_read(cov_filename)
        return cov_press

    @property
    def covariance_matrix(self):
        """Return a covariance matrix for the requested location, error if the covariance type is not available."""

        cov_filename = self.covariance_filename
        cov_press, cov_data = self._co_star_read(cov_filename)

        # "Fix" covariances that are not positive definite
        if not np.all(np.linalg.eigvals(cov_data) > 0):
            warnings.warn("Covariance matrix for species {} is not positive definite, modifying eigenvals".format(self.species))

            # Get eigen values and vector from matrix
            eigval, eigvec = np.linalg.eig(cov_data)

            # Find negative eigen values and set to the media
            eigval[np.where(eigval < 0)] = np.median(eigval)

            # Reconstruct matrix with modified eigen values
            cov_data = eigvec @ np.diag(eigval) @ np.linalg.inv(eigvec)

        return cov_data

    @property
    def constraint_pressure(self):

        con_filename = self.constraint_filename
        con_press, con_data = self._co_star_read(con_filename)
        return con_press

    @property
    def constraint_matrix(self):
        """Return a constraint matrix for the requested location, error if the covariance type is not available."""

        con_filename = self.constraint_filename
        con_press, con_data = self._co_star_read(con_filename)
        return con_data

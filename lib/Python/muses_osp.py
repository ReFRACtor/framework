import os
import re
import struct
from glob import glob
from bisect import bisect_right
import datetime as dt

import numpy as np


class MUSES_File(object):

    def __init__(self, filename, as_struct=False):

        self.filename = filename
        self.as_struct = as_struct

        self.header = {}
        self.column_names = []
        self.unit_names = []
        self.data = None

        with open(self.filename, "rb") as mf:
            self._read_file(mf)

    def _read_file(self, muses_file):

        self._read_header(muses_file)

        if 'Binary_Endian' in self.header:
            self._read_binary_data(muses_file)
        else:
            self._read_column_names(muses_file)
            self._read_ascii_data(muses_file)

    def _read_header(self, muses_file):

        # Parse heade contents until finding end of header
        curr_line = b""
        while not re.search(b"^\s*End_of_Header", curr_line):
            curr_line = muses_file.readline()
            
            header_value = re.search(b"^([\w-]+)\s*=\s*(.*)$", curr_line)
            if header_value:
                key, value = [ s.decode('UTF-8') for s in header_value.groups() ]

                # Scrub any quote marks and extra spaces
                value = value.strip('"').strip()

                # Parse datasize into an array of numbers
                data_size_re = "\s+x\s*"
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

class OSP(object):
    def __init__(self, species, base_dir, latitude, longitude, obs_time, log_cov=True, cov_dir="Covariance-scaled"):
        self.species = species
        self.base_dir = base_dir
        self.latitude = latitude
        self.longitude = longitude

        if not isinstance(obs_time, dt.datetime):
            raise Exception("obs_time must be a datetime instance")

        self.obs_time = obs_time

        self.log_cov = log_cov
        self.cov_dir = cov_dir

        self.nh3_version = "CLN"

    @property
    def climatology_base_dir(self):
        clim_base_dir = os.path.join(self.base_dir, "Climatology", self.species)
        if not os.path.exists(clim_base_dir):
            raise Exception("No climatology species directory found: {}".format(clim_base_dir))
        return clim_base_dir
 
    @property
    def pressure_filename(self):
        "Filename used for accessing the pressure grid associated with a species"

        # Find the associated pressure file
        press_file = os.path.join(self.climatology_base_dir, "Clim_Spec_{}".format(self.species))

        if not os.path.exists(press_file):
            raise Exception("Could not find pressure file: {}".format(press_file))

        return press_file

    @property
    def pressure_full_grid(self):
        "The pressure grid associated with a species converted into the format acceptable for ReFRACtor"

        press_file = self.pressure_filename
    
        # Reverse to pressure increasing order, convert hPa -> Pa
        mf = MUSES_File(press_file)
        press_data = np.flip(mf.data[:, 0] * 100)

        return press_data

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
            raise Exception("Could not determine longitude directory")

        return use_long_dir

    def _pick_lat_file(self, latitude_files):

        if len(latitude_files) == 1:
            return latitude_files[0]

        # Convert latitudes with N (North) and S (South) markings to positive and negative latitudes
        def from_ns(lat_str):
            num_val = re.sub('[NS]$', '', lat_str)
            if re.search('S', lat_str):
                return -float(num_val)
            else:
                return float(num_val)

        # Find the file that exists within the latitude range
        use_lat_file = None
        for lat_file in latitude_files:
            if re.search("NH3", lat_file):
                matches = re.search('.*(\d\d[NS])_(\d\d[NS])_{}$'.format(self.nh3_version), lat_file)
            else:
                matches = re.search('.*(\d\d[NS])_(\d\d[NS])$', lat_file)

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

        all_year_dirs = sorted([ int(y) for y in filter(lambda d: re.match('\d{4}', d), os.listdir(species_base_dir)) ])

        if len(all_year_dirs) > 0:
            year_idx = max(0, bisect_right(all_year_dirs, self.obs_time.year)-1)
            year_dir = str(all_year_dirs[year_idx])
            species_base_dir = os.path.join(species_base_dir, year_dir)
        else:
            year_dir = ''

        all_month_dirs = list(filter(lambda d: re.match('[A-Z]{3}', d), os.listdir(species_base_dir)))

        if len(all_month_dirs) == 1:
            month_dir = all_month_dirs[0]
        else:
            best_month_dir = self._month_directory_name
            if best_month_dir not in all_month_dirs:
                raise Exception("Month directory for current data {} not available at {}".format(best_month_dir, species_base_dir))

            month_dir = best_month_dir

        time_dir = "0000UT_2400UT"

        return os.path.join(year_dir, month_dir, time_dir)

    def _temporal_clim_filename(self):
        
        species_base_dir = self.climatology_base_dir

        time_rel_dir = self._time_directory_name()

        time_base_dir = os.path.join(species_base_dir, time_rel_dir)

        if not os.path.exists(time_base_dir):
            raise Exception("No climatology temporal directory found: {}".format(time_base_dir))

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
    def climatology_full_grid(self):
        "Species climatology data matching the object's time and location converted to the format handled by ReFRACtor"

        clim_file = self.climatology_filename

        if self.species == 'PSUR':
            mf = MUSES_File(clim_file, as_struct=True)

            # Extract value and convert to a matrix
            clim_data = np.array([float(mf.data['PSUR'])])

            # Convert hPa -> Pa
            clim_data *= 100
        else:
            mf = MUSES_File(clim_file)

            clim_data = mf.data[:, 0]

            # Convert to pressure increasing order
            clim_data = np.flip(clim_data)

        return clim_data

    @property
    def covariance_filename(self):
        cov_base_dir = os.path.join(self.base_dir, "Covariance", self.cov_dir)

        if not os.path.exists(cov_base_dir):
            raise Exception("Covariance base directory does not exist: {}".format(cov_base_dir))

        species_files = glob(os.path.join(cov_base_dir, "Covariance_Matrix_{}_*".format(self.species)))

        if len(species_files) == 0:
            raise Exception("No covariances files for {} found in directory: {}".format(self.species, cov_base_dir))

        if self.log_cov:
            filt_files = filter(lambda f: re.search("_Log_", f), species_files)
        else:
            filt_files = filter(lambda f: re.search("_Linear_", f), species_files)

        cov_file = self._pick_lat_file(list(filt_files))

        return cov_file

    def _covariance_read(self, cov_file):
        mf = MUSES_File(cov_file)

        species_col_names = filter(lambda cn: re.search(self.species, cn), mf.column_names)
        species_col_idx = [ mf.column_names.index(nm) for nm in species_col_names ]
        cov_data = mf.data[:, np.array(species_col_idx)]

        # Convert to pressure increasing order
        cov_data = np.flip(cov_data)

        return cov_data

    @property
    def covariance_full_grid(self):
        """Return a covariance matrix for the requested location, error if the covariance type is not available."""

        cov_filename = self.covariance_filename

        if cov_filename is None:
            raise Exception("No {} covariance file found".format(self.log_cov and "log" or "linear"))
        else:
            return self._covariance_read(cov_filename)

    def covariance_compute(self):
        """Create a simple covariance from climatology data"""

        clim_filename = self.climatology_filename

        # Replace parts of path with wildcards

        # Replace month
        clim_glob = re.sub(self._month_directory_name, "*", clim_filename)

        # Replace longitude portion
        clim_glob = re.sub("\d{3}[EW]_\d{3}[EW]", "*", clim_glob)

        # Replace latitude portion
        clim_glob = re.sub("\d{2}[NS]_\d{2}[NS]", "*", clim_glob)

        cov_inp_filenames = glob(clim_glob)
    
        inp0 = MUSES_File(cov_inp_filenames[0])

        clim_data = np.empty((inp0.data.shape[0], len(cov_inp_filenames)), dtype=float)
        clim_data[:, 0] = inp0.data[:, 0]

        for idx, inp_file in enumerate(cov_inp_filenames[1:]):
            inpD = MUSES_File(inp_file)
            clim_data[:, idx+1] = inpD.data[:, 0]

        cov_data = np.cov(clim_data)

        # Convert to pressure increasing order
        cov_data = np.flip(cov_data)

        # Add a small value to the diagonal values to make the matrix positive definite 
        # as due to numerical issues the output isn't always such
        fudge = np.mean(cov_data) * 1e-10

        diag_idx = np.diag_indices(cov_data.shape[0])
        cov_data[diag_idx] += fudge
        
        return cov_data

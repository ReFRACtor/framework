import os
import re
from glob import glob
import datetime as dt

import numpy as np

def read_muses_file(filename, as_matrix=True, **kwargs):
    "Reads a MUSES style text file"

    data = None
    with open(filename, "r") as pf:
        # Find end of header
        curr_line = ""
        while curr_line.find("End_of_Header") < 0:
            curr_line = pf.readline()
        
        # Read column names from first line after end of header
        col_names = pf.readline().strip().split()
        
        # Skip column unit labels, then read data
        data = np.genfromtxt(pf, names=col_names, skip_header=1, **kwargs)
        
    if as_matrix:
        # Convert from structured array
        mat = np.empty((data.shape[0], len(data.dtype.names)), data.dtype[0])
        for col_idx, name in enumerate(data.dtype.names):
            mat[:, col_idx] = data[name]
        return mat
    else:
        return data

class OSP(object):
    def __init__(self, base_dir, latitude, longitude, obs_time, cov_dir="Covariance-scaled"):
        self.base_dir = base_dir
        self.latitude = latitude
        self.longitude = longitude

        if not isinstance(obs_time, dt.datetime):
            raise Exception("obs_time must be a datetime instance")

        self.obs_time = obs_time

        self.cov_dir = cov_dir

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
            matches = re.search('.*(\d\d[NS])_(\d\d[NS])$', lat_file)
            beg, end = [ from_ns(lat_str) for lat_str in matches.groups() ]
            
            if self.latitude >= beg and self.latitude <= end:
                use_lat_file = lat_file
                break

        return use_lat_file

    @property
    def _month_directory_name(self):
        month_dir = self.obs_time.strftime("%b").upper()
        return month_dir

    def pressure_filename(self, species):
        "Filename used for accessing the pressure grid associated with a species"

        species_base_dir = os.path.join(self.base_dir, "Climatology", species)

        # Find the associated pressure file
        press_file = os.path.join(species_base_dir, "Clim_Spec_{}".format(species))

        if not os.path.exists(press_file):
            raise Exception("Could not find pressure file: {}".format(press_file))

        return press_file

    def pressure(self, species):
        "The pressure grid associated with a species converted into the format acceptable for ReFRACtor"

        press_file = self.pressure_filename(species)
    
        # Reverse to pressure increasing order, convert hPa -> Pa
        press_data = np.flip(read_muses_file(press_file)[:, 0] * 100)

        return press_data

    def climatology_filename(self, species):
        "Filename with the species climatology values matching the object's time and location"

        month_dir = self._month_directory_name
        time_dir = "0000UT_2400UT"

        species_base_dir = os.path.join(self.base_dir, "Climatology", species, month_dir, time_dir)

        if not os.path.exists(species_base_dir):
            raise Exception("No climatology directory found: {}".format(species_base_dir))

        long_dir = self._pick_long_dir(os.listdir(species_base_dir))
        lat_file = self._pick_lat_file(os.listdir(os.path.join(species_base_dir, long_dir)))

        clim_file = os.path.join(species_base_dir, long_dir, lat_file)

        return clim_file

    def climatology(self, species):
        "Species climatology data matching the object's time and location converted to the format handled by ReFRACtor"

        clim_file = self.climatology_filename(species)

        if species == 'PSUR':
            clim_data = read_muses_file(clim_file, as_matrix=False)

            # Extract value and convert to a matrix
            clim_data = np.array([float(clim_data['PSUR'])])

            # Convert hPa -> Pa
            clim_data *= 100
        else:
            clim_data = read_muses_file(clim_file)[:, 0]

            # Convert to pressure increasing order
            clim_data = np.flip(clim_data)

        return clim_data

    def covariance_filename(self, species):
        cov_base_dir = os.path.join(self.base_dir, "Covariance", self.cov_dir)

        if not os.path.exists(cov_base_dir):
            raise Exception("Covariance base directory does not exist: {}".format(cov_base_dir))

        species_files = glob(os.path.join(cov_base_dir, "Covariance_Matrix_{}_*".format(species)))

        if len(species_files) == 0:
            raise Exception("No covariances files for {} found in directory: {}".format(species, cov_base_dir))

        cov_file = self._pick_lat_file(species_files)

        return cov_file

    def _covariance_read(self, species):
        cov_file = self.covariance_filename(species)

        cov_data = read_muses_file(cov_file)[:, 2:]

        # Convert to pressure increasing order
        cov_data = np.flip(cov_data)

        return cov_data

    def _covariance_compute(self, species):

        clim_filename = self.climatology_filename(species)

        # Replace parts of path with wildcards

        # Replace month
        clim_glob = re.sub(self._month_directory_name, "*", clim_filename)

        # Replace longitude portion
        #clim_glob = re.sub("\d{3}[EW]_\d{3}[EW]", "*", clim_glob)

        # Replace latitude portion
        #clim_glob = re.sub("\d{2}[NS]_\d{2}[NS]", "*", clim_glob)

        print(clim_glob)

        cov_inp_filenames = glob(clim_glob)
    
        inp0 = read_muses_file(cov_inp_filenames[0]) 

        clim_data = np.empty((inp0.shape[0], len(cov_inp_filenames)), dtype=float)
        clim_data[:, 0] = inp0[:, 0]

        for idx, inp_file in enumerate(cov_inp_filenames[1:]):
            inpD = read_muses_file(inp_file)
            clim_data[:, idx+1] = inpD[:, 0]

        cov_data = np.cov(clim_data)

        # Convert to pressure increasing order
        cov_data = np.flip(cov_data)

        # Add a small value to the diagonal values to make the matrix positive definite 
        # as due to numerical issues the output isn't always such
        fudge = np.mean(cov_data) * 1e-10

        diag_idx = np.diag_indices(cov_data.shape[0])
        cov_data[diag_idx] += fudge
        
        return cov_data

    def covariance(self, species, in_log=False, force_compute=False):
        """Return a covariance matrix for the requested location, returns either a linear or log matrix. 
           If one or the other is not available then a covariance will be created from climatology values"""

        cov_filename = self.covariance_filename(species)

        if (re.search('_Log_', cov_filename) and not in_log) or force_compute:
            return self._covariance_compute(species)
        else:
            return self._covariance_read(species)

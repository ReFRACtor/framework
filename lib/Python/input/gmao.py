import os
import logging

from glob import glob

import numpy as np
import netCDF4

logger = logging.getLogger(__name__)

class GMAOFinder(object):
    "Find a GMAO filename from a directory structure based on time"

    search_time = '{when.year:04d}{when.month:02d}{when.day:02d}_{hour:02d}00'
    search_version = 'V01'
    search_extension = '.nc4'

    # Find the first filename matching these patterns
    search_names_3d = [ 'asm.inst3_3d_asm_Np', 'fcst.inst3_3d_asm_Np' ]
    search_names_2d = [ 'asm.inst3_2d_asm_Nx', 'fcst.tavg1_2d_slv_Nx' ]

    # Directory organizations to search
    search_directories = [ '{base_dir}',
                           '{base_dir}/{when.year:04d}',
                           '{base_dir}/{when.year:04d}/{when.month:02d}', 
                           '{base_dir}/{when.year:04d}/{when.month:02d}/{when.day:02d}', ]

    def __init__(self, when, base_dir, hour_delta=3, search_directories=None, **kwargs):

        self.when = when
        self.base_dir = base_dir.rstrip("/")

        if search_directories is not None:
            self.search_directories = search_directories

        # Filenames are every N hours
        self.hour = int(np.round(self.when.hour / hour_delta, 0) * hour_delta)

    def _find_file(self, search_names):

        for srch_name in search_names:
            filename_glob = f"*.{srch_name}.*.{self.search_time}*.{self.search_version}{self.search_extension}".format(when=self.when, hour=self.hour)
            
            for curr_dir in self.search_directories:
                directory_glob = curr_dir.format(base_dir=self.base_dir, when=self.when)

                search_glob = os.path.join(directory_glob, filename_glob)

                logger.debug(f"Searching for GMAO using pattern: {search_glob}")

                filenames_found = glob(search_glob)

                if len(filenames_found) == 1:
                    logger.debug(f"Found GMAO file for time {self.when}: {filenames_found[0]}")
                    return filenames_found[0]
                elif len(filenames_found) > 1:
                    raise Exception(f"Too many GMAO filenames found ({len(filenames_found)}) using glob: {search_glob}")
        
        raise Exception(f"No filenames found from base path {self.base_dir} for file name types: {search_names}")

    @property
    def filename_3d(self):

        return self._find_file(self.search_names_3d)

    @property
    def filename_2d(self):

        return self._find_file(self.search_names_2d)

class GMAOReader(object):
    "Directly read the indicated GMAO files"

    def __init__(self, longitude, latitude, filename_2d, filename_3d, **kwargs):

        self.longitude = longitude
        self.latitude = latitude

        self.filename_2d = filename_2d
        self.filename_3d = filename_3d

        # Nominally these should be the same, but lets not make assumptions
        self.idx_lon_2d, self.idx_lat_2d = self._closest_index(self.filename_2d)
        self.idx_lon_3d, self.idx_lat_3d = self._closest_index(self.filename_3d)

        if self.idx_lon_2d != self.idx_lon_3d:
            logger.warning(f"GMAO 2D file {filename_2d} and 3D file {filename_3d} have different longitude indexes")

        if self.idx_lat_2d != self.idx_lat_3d:
            logger.warning(f"GMAO 2D file {filename_2d} and 3D file {filename_3d} have different latitude indexes")

    def _closest_index(self, filename):

        with netCDF4.Dataset(filename, "r") as gmao_data:

            gmao_lon = gmao_data['lon'][:]
            gmao_lat = gmao_data['lat'][:]
            
            lon_diff = np.abs(gmao_lon - self.longitude)
            lon_idx = np.where(lon_diff == np.min(lon_diff))[0][0]

            lat_diff = np.abs(gmao_lat - self.latitude)
            lat_idx = np.where(lat_diff == np.min(lat_diff))[0][0]

        return lon_idx, lat_idx

    def _get_2d_value(self, name):

        with netCDF4.Dataset(self.filename_2d) as gmao_2d:

            value = gmao_2d[name][..., self.idx_lat_2d, self.idx_lon_2d]
            
            # Return only valid values and convert to normal numpy array
            return np.array(value[value.mask == False])

    def _get_3d_value(self, name):

        with netCDF4.Dataset(self.filename_3d) as gmao_3d:

            value = gmao_3d[name][..., self.idx_lat_3d, self.idx_lon_3d]

            # Return only valid values and convert to normal numpy array
            return np.array(value[value.mask == False])

    @property
    def pressure(self):

        # Truncate pressure to the same grid as used by the humidity
        with netCDF4.Dataset(self.filename_3d) as gmao_3d:
            mask = gmao_3d['QV'][0, :, self.idx_lat_3d, self.idx_lon_3d].mask

            # Put in pressure increasing order and convert to Pascals
            return np.array(gmao_3d['lev'][mask == False][::-1]) * 100

    @property
    def temperature(self):

        # Put in pressure increasing order
        return self._get_3d_value("T")[::-1] 

    @property
    def h2o(self):

        # convert specified humidity to vmr
        Md=28.966 # molecular mass of air
        Mw=18.016 # molecular mass of water
        qv_2_vmr = Md/Mw

        # Put in pressure increasing order
        return self._get_3d_value("T")[::-1] * qv_2_vmr

    @property
    def surface_pressure(self):

        return self._get_2d_value("SLP").item()

    @property
    def surface_temperature(self):

        return self._get_2d_value("TS").item()

    @property
    def tropopause_pressure(self):

        return self._get_2d_value("TROPPT").item()

class GMAOSupplier(object):
    """Loads GMAO values given latitude, longitude and time. Ensures that the pressure grid includes
    the surface pressure"""

    def __init__(self, longitude, latitude, when, base_dir, **kwargs):

        search_gmao = GMAOFinder(when, base_dir, **kwargs)

        self.reader = GMAOReader(longitude, latitude, search_gmao.filename_2d, search_gmao.filename_3d, **kwargs)

        press = self.reader.pressure
        psurf = self.reader.surface_pressure

        # Ensure that the grid includes surface pressure
        self.extend_levels = False
        if psurf > np.max(press):
            self.extend_levels = True

    @property
    def pressure(self):

        press = self.reader.pressure

        if self.extend_levels:        
            psurf = self.reader.surface_pressure
            return np.append(press, psurf)
        else:
            return press

    @property
    def temperature(self):

        temp = self.reader.temperature 

        if self.extend_levels:
            return np.append(temp, temp[-1])
        else:
            return temp

    @property
    def h2o(self):

        wv = self.reader.h2o

        if self.extend_levels:
            return np.append(wv, wv[-1])
        else:
            return wv

    @property
    def surface_pressure(self):

        return self.reader.surface_pressure

    @property
    def surface_temperature(self):

        return self.reader.surface_temperature

    @property
    def tropopause_pressure(self):

        return self.reader.tropopause_pressure

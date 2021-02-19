#!/usr/bin/env python3

import os
import sys
import logging
import shutil
import logging

import progressbar

import numpy as np
import netCDF4

from refractor import framework as rf
from refractor.executor import StrategyExecutor

logger = logging.getLogger(__file__)

template_basename = "optical_properties_template.nc"
if "REFRACTOR_INPUTS" in os.environ and os.path.exists(os.path.join(os.environ["REFRACTOR_INPUTS"], template_basename)):
    # From installed path
    template_filename = os.path.join(os.environ["REFRACTOR_INPUTS"], template_basename)
else:
    # From source path
    template_dir = os.path.join(os.path.dirname(__file__), '../input/template')
    template_filename = os.path.join(template_dir, template_basename)

def write_scene_optical_depth(step_index, config_inst, output, num_moment=16, num_scatt=6):

    fm = config_inst.forward_model
    atm = config_inst.atmosphere

    opt_prop_grp = output['OpticalProperties']
   
    for gas_index in range(atm.absorber.number_species):
        gas_name = atm.absorber.gas_name(gas_index)
        opt_prop_grp['gas_name'][step_index, gas_index, :len(gas_name)] = \
            netCDF4.stringtochar(np.array(gas_name, dtype=f'S{len(gas_name)}'))

    n_sensor = fm.num_channels
        
    n_grid_max = 0
    for sensor_idx in range(n_sensor):
        n_grid_max += fm.spectral_grid.high_resolution_grid(sensor_idx).data.shape[0]

    n_gas = atm.absorber.number_species
    if atm.aerosol is not None:
        n_aerosol = atm.aerosol.number_particle
    else:
        n_aerosol = 0
    n_layer = atm.number_layer

    for sensor_idx in range(n_sensor):
        logger.debug(f"Create optical properties for scenario {step_index+1}, step {sensor_idx+1}")

        grid = fm.spectral_grid.high_resolution_grid(sensor_idx)
        n_grid = grid.data.shape[0]

        opt_prop_grp['grid'][step_index, sensor_idx, :n_grid] = grid.data
        opt_prop_grp['grid_size'][step_index, sensor_idx] = n_grid

        # Store in temporary arrays first for speed, writing to file indexes individually is slow
        gas_od = np.empty((n_grid, n_gas, n_layer))
        aer_ext = np.empty((n_grid, n_aerosol, n_layer))
        ray_od = np.empty((n_grid, n_layer))
        tot_od = np.empty((n_grid, n_layer))
        tot_ssa = np.empty((n_grid, n_layer))
        tot_pf = np.empty((n_grid, num_moment, n_layer, num_scatt))

        bar = progressbar.ProgressBar(max_value=n_grid)
            
        for grid_idx, grid_value in enumerate(grid.data):
            grid_point = rf.DoubleWithUnit(grid_value, grid.units)
            opt_prop = atm.optical_properties(grid_point.convert_wave("cm^-1").value, 0)
            
            gas_od[grid_idx, ...] = opt_prop.gas_optical_depth_per_particle().value.transpose()
            ray_od[grid_idx, ...] = opt_prop.rayleigh_optical_depth().value
            tot_od[grid_idx, ...] = opt_prop.total_optical_depth().value
            tot_ssa[grid_idx, ...] = opt_prop.total_single_scattering_albedo().value

            if atm.aerosol is not None:
                aer_ext[grid_idx, ...] = opt_prop.aerosol_extinction_optical_depth_per_particle().value.transpose()

                pf = opt_prop.total_phase_function_moments().value

                tot_pf[grid_idx, ...] = pf[:num_moment, :, :num_scatt]
     
            bar.update(grid_idx+1)

        bar.finish()

        logger.debug("Saving standard optical properties to file")

        opt_prop_grp['gas_optical_depth_per_particle'][step_index, sensor_idx, :n_grid, :n_gas, :n_layer] = gas_od
        opt_prop_grp['rayleigh_optical_depth'][step_index, sensor_idx, :n_grid, :n_layer] = ray_od
        opt_prop_grp['total_optical_depth'][step_index, sensor_idx, :n_grid, :n_layer] = tot_od
        opt_prop_grp['total_single_scattering_albedo'][step_index, sensor_idx, :n_grid, :n_layer] = tot_ssa

        if atm.aerosol is not None:
            logger.debug("Saving aerosol optical properties to file")
            opt_prop_grp['aerosol_extinction_optical_depth_per_particle'][step_index, sensor_idx, :n_grid, :n_aerosol, :n_layer] = aer_ext
            opt_prop_grp['total_phase_function_moments'][step_index, sensor_idx, :n_grid, :num_moment, :n_layer, :num_scatt] = tot_pf

def write_optical_properties(exc, output_filename, step_indexes=None, **kwargs):

    if step_indexes is None:
        step_indexes = [0]

    ## Write Output

    # Copy template to output
    logger.debug(f"Creating optical properties output file: {output_filename} from template file: {template_filename}")
    shutil.copyfile(template_filename, output_filename)
    output = netCDF4.Dataset(output_filename, "a")

    for step_index in step_indexes:
        step_keywords = exc.strategy_list[step_index]
        config_inst = exc.config_instance(**step_keywords)

        if 'step_id' in config_inst.file_config:
            id_len = len(config_inst.file_config.step_id)
            output['Scenario/scene_id'][step_index, :id_len] = netCDF4.stringtochar(np.array(config_inst.file_config, dtype=f'S{id_len}'))

        write_scene_optical_depth(step_index, config_inst, output, **kwargs)

    output.close()

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Write optical properties set up by configuration")
    
    parser.add_argument("config_filename", metavar="FILE",
        help="File containing a configuration method returning a configuration instance")

    parser.add_argument("--strategy_filename", "-s", metavar="FILE", default=None,
        help="File containing the strategy defining steps to take for the retrival, if not defined a single step will be taken")

    parser.add_argument("--output_filename", "-o", metavar="FILE", default=None,
        help="File to write results of configuration execution")

    parser.add_argument("--step_index", metavar="INT", type=int, action="append", default=None,
        help="Strategy step indexes to us instead of all. Specify multiple times for multiple")

    parser.add_argument("--num_moment", metavar="INT", type=int, default=16,
        help="Number of aerosol moments to save from the phase function")

    parser.add_argument("--verbose", "-v", action="store_true", default=False,
        help="Increase level of verbosity")

    args = parser.parse_args()
 
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stdout)

    exc = StrategyExecutor(args.config_filename, 
                           strategy_filename=args.strategy_filename)

    write_optical_properties(exc, args.output_filename, step_indexes=args.step_index, num_moment=args.num_moment)

if __name__ == "__main__":
    main()

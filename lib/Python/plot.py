import numpy as np
import matplotlib.pyplot as plt

def plot_absorbers(atm):
    "Plot all absorbers from an atmosphere class"

    alt_grid = atm.altitude(0).value.value
    for sidx in range(atm.absorber.number_species):
        gas_name = atm.absorber.gas_name(sidx)
        avmr = atm.absorber.absorber_vmr(gas_name)
        
        # Get the vmr data on the same grid as the altitude grid
        vmr_profile = avmr.vmr_grid(atm.pressure).value
        
        plt.figure()
        plt.plot(vmr_profile, alt_grid)
        plt.title(gas_name)
        plt.ylabel("Altitude (km)")
        plt.xlabel("VMR (ppm)")

def mw_ranges(domain, grid_jump_tol=1):
    "Separate consecutive microwindows into pieces"

    dom_diff = (domain[1:] - domain[:-1])
    wend = np.where(dom_diff > grid_jump_tol)[0]
    wbeg = np.concatenate(([0], wend+1))
    wend = np.concatenate((wend, [domain.shape[0]-1]))

    return zip(wbeg, wend)

def plot_microwindows(meas_grid, meas_rad, mod_grid, mod_rad, grid_units=None, rad_units=None, meas_uncert=None, show_residual=False, title_prefix=""):

    for beg_idx, end_idx in mw_ranges(meas_grid):
        plt.figure()
        if show_residual:
            plt.plot(meas_grid[beg_idx:end_idx+1], 
                (mod_rad[beg_idx:end_idx+1] - meas_rad[beg_idx:end_idx+1]) / uncert[beg_idx:end_idx+1] )
        else:
            plt.plot(meas_grid[beg_idx:end_idx+1], meas_rad[beg_idx:end_idx+1])
            plt.plot(mod_grid[beg_idx:end_idx+1], mod_rad[beg_idx:end_idx+1])

            if rad_units is not None:
                plt.ylabel("$ {} $".format(rad_units))
            plt.legend(['measured', 'modeled'], loc='best')

        if grid_units is not None:
            plt.xlabel("$ {} $".format(grid_units))
        plt.title("{}MW {} - {}, {} - {}".format(title_prefix, meas_grid[beg_idx], meas_grid[end_idx], beg_idx, end_idx))
    
    if meas_uncert is not None:
        plt.figure()
        plt.plot((mod_rad - meas_rad)/uncert)
        plt.title("{}Residual All Samples".format(title_prefix))

        plt.figure()
        plt.plot(meas_rad/meas_uncert)
        plt.title("{}~ SNR".format(title_prefix))

        plt.figure()
        plt.plot(uncert)
        plt.title("{}Uncertainty All Samples".format(title_prefix))

def plot_fm(fm, l1b, channel_idx=0, show_residual=False, rad=None):
    pix_list = np.array(fm.spectral_grid.pixel_list(channel_idx))

    if rad is None:
        rad = fm.radiance(channel_idx)
    
    l1b_grid = l1b.sample_grid(channel_idx).data[pix_list]
    l1b_rad = l1b.radiance(channel_idx).data[pix_list]
    uncert = l1b.radiance(channel_idx).uncertainty[pix_list]

    grid_units = l1b.sample_grid(channel_idx).units.name
    meas_units = l1b.radiance(channel_idx).units.name
    
    mod_grid = rad.spectral_domain.data.copy()
    mod_rad = rad.spectral_range.data.copy()

    plot_microwindows(l1b_grid, l1b_rad, mod_grid, mod_rad, grid_units, meas_units, uncert, show_residual=show_residual)

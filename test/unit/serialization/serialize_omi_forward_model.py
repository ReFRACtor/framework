#!/usr/bin/env python

# This create a serialization of an example omi forward model. We extracted
# this from a real test run 20160414_23_394_23. This is really just meant
# to be "an example", so the specific sounding isn't all that important.

from refractor import framework as rf
from refractor import write_shelve

from refractor.input.paths import cross_section_filenames, cross_section_file_conversion_factors

import numpy as np
import h5py
import os

input_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "../data/in")) + "/"
serialized_fname = input_dir + "forward_model/omi_underlying_forward_model.xml"

# We extracted omi data, and hardcode it here so we don't depend on omi
# to regenerate this serialized model.
#
# These values were just read using the code:
#
# from omi import *
# from refractor import framework as rf
# import numpy as np
# t = RefractorUip("/home/smyth/output_py/omi/2016-04-14/setup-targets/Global_Survey/20160414_23_394_23/Table.asc")
# mp = UipToRefactorMap(t.uip)
# c = OmiObjectCreator(t.uip,0)
# np.set_printoptions(precision=20)

# mp.across_track_indexes(0)[0]
across_track_index = 23

# mp.earth_sun_distance
earth_sun_distance = 150089039872.0

# mp.micro_windows(0).value
micro_windows = np.array([[346.5, 347.5]])

# c.instrument_spectral_domain.data
instrument_spectral_domain_data = np.array(
    [346.59747645990745, 346.7403556701846 , 346.88317556485816,
     347.02593618744083, 347.16863758104216, 347.31127978836787,
     347.45386285172077])

# mp.atmosphere_column("pressure") * 100
plev = np.array(
    [1.0000000000000000e+01, 2.1544300000000000e+01,
     3.8311799999999998e+01, 6.8129199999999997e+01,
     1.0000000000000000e+02, 1.3335200000000000e+02,
     1.6156000000000000e+02, 2.1544300000000001e+02,
     2.6101600000000002e+02, 3.1622699999999998e+02,
     4.6416000000000003e+02, 5.1089799999999997e+02,
     6.8129099999999994e+02, 8.2540200000000004e+02,
     9.0851400000000012e+02, 1.0000000000000000e+03,
     1.1007000000000000e+03, 1.2115300000000000e+03,
     1.3335200000000000e+03, 1.4678000000000002e+03,
     1.6155999999999999e+03, 1.7782800000000002e+03,
     1.9573399999999999e+03, 2.1544299999999998e+03,
     2.3713600000000001e+03, 2.6101700000000001e+03,
     2.8729900000000002e+03, 3.1622900000000000e+03,
     3.4807100000000000e+03, 3.8311900000000001e+03,
     4.2169600000000000e+03, 4.6415799999999999e+03,
     5.1089600000000000e+03, 5.6233899999999994e+03,
     6.1896299999999992e+03, 6.8129499999999989e+03,
     7.4989599999999991e+03, 8.2540599999999995e+03,
     9.0851800000000003e+03, 1.0000000000000000e+04,
     1.1006900000000000e+04, 1.2115200000000001e+04,
     1.3335200000000001e+04, 1.4677900000000000e+04,
     1.6156100000000000e+04, 1.7782900000000001e+04,
     1.9573500000000000e+04, 2.1544399999999998e+04,
     2.3713700000000001e+04, 2.6101600000000002e+04,
     2.8729799999999999e+04, 3.1622699999999997e+04,
     3.4806900000000001e+04, 3.8311700000000004e+04,
     4.2169799999999996e+04, 4.6416000000000000e+04,
     5.1089800000000003e+04, 5.6234199999999997e+04,
     6.1896599999999999e+04, 6.8129100000000006e+04,
     7.4989300000000003e+04, 8.2540200000000012e+04,
     9.0851399999999994e+04, 1.0000000000000000e+05,
     1.0147354687500003e+05])

# mp.omi_params["cloud_pressure"] * 100.0
cloud_pressure = 65500.0

# mp.atmosphere_column("TATM")
temperature_levels = np.array(
    [227.09063720703125, 238.06014518011747, 247.9513508952482 ,
     255.45985665566292, 252.05914306640625, 247.81659276701745,
     244.9881821673292 , 241.44002596182148, 240.4029340659358 ,
     238.91132912433085, 232.9867434126154 , 231.29077037580797,
     226.18198665229613, 224.41776735108155, 223.67060445641607,
     222.92340087890625, 223.0570635343525 , 223.19071412774932,
     223.32436521762182, 223.45802268716176, 223.5916789100213 ,
     223.72533301080605, 223.8589865589824 , 223.07685269163323,
     222.02927852743204, 220.9815484320261 , 219.9339687734463 ,
     219.05707593300883, 218.32040612970977, 217.5837415734727 ,
     217.03191042245126, 216.63095888131625, 216.20113885493686,
     215.6717339085364 , 215.1423190168783 , 214.61285294111053,
     214.57427054214318, 214.7287647530646 , 214.8832564805342 ,
     215.03775024414062, 215.68044262838416, 216.3231434390266 ,
     216.96589636412992, 217.60857903239634, 216.28642561078664,
     214.39008748263186, 212.4937561712388 , 213.1970385598021 ,
     214.65362920333132, 218.19508337616114, 224.29002583967434,
     230.53981263587198, 236.91667086084584, 242.25954916368235,
     246.5540877988173 , 250.1889395480286 , 254.5180214696311 ,
     260.19581461861503, 265.5544780044012 , 269.78800660563746,
     273.14080883316444, 272.72152807894287, 279.0609940995202 ,
     286.82574462890625, 286.82574462890625])

# mp.atmosphere_column("O3")
o3_vmr_levels = np.array(
    [7.2132046260759609e-07, 1.0912349576548504e-06,
     1.4885425799767671e-06, 2.0305073336902930e-06,
     3.0388817017368270e-06, 4.1119249342401622e-06,
     5.0303798445834408e-06, 5.7386008770588536e-06,
     6.2653299262418214e-06, 6.8403846940672217e-06,
     8.1537243750512442e-06, 8.0142765809954363e-06,
     7.6100832447428538e-06, 7.5420257399991535e-06,
     7.5082258354330517e-06, 7.4172656409340817e-06,
     7.3274061278404856e-06, 7.2386432297123408e-06,
     7.0922125886425295e-06, 6.9487373154825623e-06,
     6.8081658246209090e-06, 6.5102670037305346e-06,
     6.2254042205718432e-06, 5.9530059922994415e-06,
     5.4366554982230307e-06, 4.9650250727046469e-06,
     4.5343677485424209e-06, 4.0309301905529154e-06,
     3.5834023071197631e-06, 3.1855632736432420e-06,
     2.7638746303455839e-06, 2.3980023546645406e-06,
     2.0805608342232576e-06, 1.6764810216750137e-06,
     1.3508748093858100e-06, 1.0884851725418538e-06,
     8.7114396521552429e-07, 6.9719717340910282e-07,
     5.5798542001817292e-07, 4.3997029512579159e-07,
     3.4691848441742429e-07, 2.7354586831260635e-07,
     2.1798788353631903e-07, 1.7371819360112942e-07,
     1.3843274907148251e-07, 1.1576445135856372e-07,
     9.6808138673089570e-08, 8.0955901863866402e-08,
     6.9906741888590096e-08, 6.0364998463238286e-08,
     5.4051661316549779e-08, 4.8398499029885246e-08,
     4.5052733532417181e-08, 4.1938282152302465e-08,
     4.0050832153006740e-08, 3.8248494338833188e-08,
     3.6807463767747135e-08, 3.5420730559872555e-08,
     3.3879917115943875e-08, 3.2406143996314936e-08,
     2.9752521396487041e-08, 2.7316218334316160e-08,
     2.5164498122494275e-08, 2.3182166907932117e-08,
     2.2893973089497575e-08])

# mp.muses_fm_spectral_domain(0).data
muses_fm_spectral_domain_data = np.array(
    [343.58322979298174, 343.72736539665163, 343.87144067778354,
     344.0154556883598 , 344.1594104799597 , 344.30330510375916,
     344.4471396105307 , 344.5909140506435 , 344.7346284740635 ,
     344.87828293035335, 345.02187746867213, 345.1654121377759 ,
     345.30888698601706, 345.4523020613451 , 345.59565741130587,
     345.7389530830419 , 345.8821891232926 , 346.02536557839386,
     346.16848249427835, 346.31153991647534, 346.4545378901108 ,
     346.59747645990745, 346.7403556701846 , 346.88317556485816,
     347.02593618744083, 347.16863758104216, 347.31127978836787,
     347.45386285172077, 347.59638681300027, 347.7388517137024 ,
     347.8812575949199 , 348.02360449734203, 348.16589246125505,
     348.3081215265415 , 348.45029173268097, 348.5924031187494 ,
     348.73445572341956, 348.87644958496094, 349.01838474123963,
     349.16026122971834, 349.3020790874566 , 349.44383835111057,
     349.5855390569329 , 349.7271812407731 , 349.86876493807745,
     350.0102901838886 , 350.15175701284625, 350.2931654591863 ,
     350.43451555674176])

# mp.latitude_with_unit(0)
# mp.surface_height_with_unit(0)
# float(mp.solar_zenith(0))
# float(mp.observation_zenith(0))
# float(mp.relative_azimuth(0))
# mp.omi_params['surface_albedo_uv2']
# mp.omi_params['surface_albedo_slope_uv2']
# mp.omi_params['cloud_Surface_Albedo']
# mp.omi_params['ring_sf_uv2']

latitude_with_unit = rf.DoubleWithUnit(-36.95997619628906, "deg")
surface_height_with_unit = rf.DoubleWithUnit(0, "m")
sza = np.array([54.7393913269043])
oza = np.array([15.271627426147461])
raz = np.array([62.39064407348633])
sza_with_unit = rf.ArrayWithUnit(sza, "deg")
oza_with_unit = rf.ArrayWithUnit(oza, "deg")
raz_with_unit = rf.ArrayWithUnit(raz, "deg")
surface_albedo = 0.088
surface_albedo_slope = 0
cloud_albedo = 0.8
raman_scale_factor = 1.9

# --------------
# Solar model
solar_reference_fname = input_dir + "solar/omisol_v003_avg_nshi_backup.h5"
f = h5py.File(solar_reference_fname, "r")
filter_name = "UV2"
wav_vals = f[f"WAV_{filter_name}"][:, across_track_index]
irad_vals = f[f"SOL_{filter_name}"][:, across_track_index]
one_au = 149597870691
irad_vals *= (one_au / earth_sun_distance) ** 2
# File does not have units contained within it
# Same units as the OMI L1B files, but use irradiance units
sol_domain = rf.SpectralDomain(wav_vals, rf.Unit("nm"))
sol_range = rf.SpectralRange(irad_vals, rf.Unit("ph / nm / s"))
sol_spec = rf.Spectrum(sol_domain, sol_range)
solar_model = rf.SolarReferenceSpectrum(sol_spec, None)

# --------------
# spec_win
spec_win = rf.SpectralWindowRange(rf.ArrayWithUnit(np.array([micro_windows]), "nm"))

# --------------
# Spectrum sampling
spectrum_sampling = rf.IdentitySpectrumSampling(1)

# --------------
# Instrument correction
instrument_correction = rf.vector_vector_instrument_correction()
instrument_correction.push_back(rf.vector_instrument_correction())

# --------------
# instrument spectral domain
instrument_spectral_domain = rf.SpectralDomain(instrument_spectral_domain_data,
                                               rf.Unit("nm"))

# --------------
# Instrument

sg = rf.SampleGridSpectralDomain(instrument_spectral_domain, "UV2")
ils_vec = rf.vector_ils()
ils_vec.push_back(rf.IdentityIls(sg))
instrument = rf.IlsInstrument(ils_vec, instrument_correction)

# --------------
# Pressure
pressure_clear = rf.PressureSigma(plev, plev[-1])
pressure = rf.PressureWithCloudHandling(pressure_clear, cloud_pressure)

# --------------
# Temperature
temperature = rf.TemperatureLevel(temperature_levels, pressure_clear)

# --------------
# altitude
altitude = rf.vector_altitude()
altitude.push_back(rf.AltitudeHydrostatic(pressure, temperature,
                                          latitude_with_unit,
                                          surface_height_with_unit))

# --------------
# Rayleigh model
constants = rf.DefaultConstant()
rayleigh = rf.RayleighBodhaine(pressure, altitude, constants)

# --------------
# Absorber
vmrs = rf.vector_absorber_vmr()
vmrs.push_back(rf.AbsorberVmrLevel(pressure_clear, o3_vmr_levels,
                                   "O3", rf.StateMappingLog()))

xsec_data = np.loadtxt(cross_section_filenames["O3"])
spec_grid = rf.ArrayWithUnit(xsec_data[:, 0], "nm")
xsec_values = rf.ArrayWithUnit(xsec_data[:, 1:], "cm^2")

o3_table = rf.XSecTableTempDep(spec_grid, xsec_values, cross_section_file_conversion_factors["O3"])

xsec_tables = rf.vector_xsec_table()
xsec_tables.push_back(o3_table)

absorber = rf.AbsorberXSec(vmrs, pressure, temperature, altitude, xsec_tables)

# --------------
# Ground

albedo = np.zeros((1, 2))
band_reference = rf.ArrayWithUnit(np.array([(386+306)/ 2,]),"nm")
albedo[0, 0] = surface_albedo
which_retrieved = np.array([[True, True],], dtype = bool)
ground_clear = rf.GroundLambertian(albedo, band_reference, ["UV2"],
               rf.StateMappingAtIndexes(np.ravel(which_retrieved)))
ground = rf.GroundWithCloudHandling(ground_clear, cloud_albedo, False)

# --------------
# Relative humidity
relative_humidity = rf.RelativeHumidity(absorber, temperature, pressure)

# --------------
# Atmosphere
atmosphere = rf.AtmosphereStandard(absorber, pressure, temperature, rayleigh,
                   relative_humidity, ground, altitude, constants)

# --------------
# RamanSiorisEffect
muses_fm_spectral_domain = rf.SpectralDomain(muses_fm_spectral_domain_data,
                                             rf.Unit("nm"))
raman_effect = rf.RamanSiorisEffect(muses_fm_spectral_domain,
           raman_scale_factor, 0, sza_with_unit[0], oza_with_unit[0],
           raz_with_unit[0], atmosphere, solar_model, rf.StateMappingLinear())

# --------------
# Spectrum effect
spectrum_effect = rf.vector_vector_spectrum_effect()
per_channel_eff = rf.vector_spectrum_effect()
per_channel_eff.push_back(raman_effect)
spectrum_effect.push_back(per_channel_eff)

# --------------
# Radiative transfer
num_streams = 4
num_mom = 16
use_thermal_emission = False
use_solar_sources = True
pure_nadir = False
multiple_scattering_only = False
use_thermal_scattering = True
a = np.zeros((1, 4))
a[:,0] = 1
stokes = rf.StokesCoefficientConstant(a)
radiative_transfer = rf.LidortRt(atmosphere, stokes, sza, oza, raz,
                                 pure_nadir, num_streams, num_mom,
                                 multiple_scattering_only, use_solar_sources,
                                 use_thermal_emission, use_thermal_scattering)
# Change RT flags to match py-retrieve
lid_interface = radiative_transfer.rt_driver.lidort_interface
lid_interface.lidort_modin.mbool().ts_do_sscorr_nadir(False)
lid_interface.lidort_modin.mbool().ts_do_sscorr_outgoing(False)
lid_interface.lidort_modin.mbool().ts_do_rayleigh_only(True)
lid_interface.lidort_modin.mbool().ts_do_double_convtest(False)
lid_interface.lidort_modin.mbool().ts_do_deltam_scaling(False)
lid_interface.lidort_modin.mchapman().ts_earth_radius(6371.0)
lid_interface.lidort_fixin.cont().ts_nstreams(2)
lid_interface.lidort_modin.mcont().ts_nmoments_input(2)

# --------------
# Underlying forward model (before cloud fraction)
underlying_forward_model = rf.StandardForwardModel(instrument,
       spec_win, radiative_transfer, spectrum_sampling, spectrum_effect)
underlying_forward_model.setup_grid()

write_shelve(serialized_fname, underlying_forward_model)

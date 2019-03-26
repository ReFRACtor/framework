require "example_base_config"

config = ExampleBaseConfig:new()
config.fm.atmosphere.absorber.gases = {"CO2", "H2O", "O2", "CF4"}

-- This has 2 broadeners
config.fm.atmosphere.absorber.O2.absco_aer = "absco_aer/O2_01205-01755_v0.0_init.nc"
config.fm.atmosphere.absorber.O2.absco_aer_interpolation_type=1

-- This has no broadeners
config.fm.atmosphere.absorber.CF4 = {
   absco_aer = "absco_aer/CF4_01205-01755_v0.0_init.nc",
   table_scale = 1.0,
   creator = ConfigCommon.vmr_level_constant_well_mixed,
   absco_aer_interpolation_type=1,
   -- Nonsense value, we just need something so we copy the O2 apriori
   apriori = ConfigCommon.hdf_read_double_1d("Gas/O2/average_mole_fraction")
}



config:do_config()


require "example_base_config"

config = ExampleBaseConfig:new()
config.fm.atmosphere.absorber.O2.absco_aer = "absco_aer/O2_06140-13230_v0.0_init.nc"
config.fm.atmosphere.absorber.O2.absco_aer_interpolation_type=1

config:do_config()


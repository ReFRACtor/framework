------------------------------------------------------------
--- Do a run using the ExampleBaseConfig, but with BreonVeg ground
------------------------------------------------------------

require "example_base_config"

config = ExampleBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_brdf_veg

config:do_config()

------------------------------------------------------------
--- Do a run using the ExampleBaseConfig, but with BreonSoil ground
------------------------------------------------------------

require "example_base_config"

config = ExampleBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_brdf_soil

config:do_config()

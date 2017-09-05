------------------------------------------------------------
--- Do a run using the ExampleBaseConfig, but with Coxmunk ground
------------------------------------------------------------

require "example_base_config"

config = ExampleBaseConfig:new()

config.fm.atmosphere.ground.creator = ConfigCommon.ground_coxmunk

config:do_config()

------------------------------------------------------------
--- Do a run using the ExampleBaseConfig, but with Coxmunk+Lamb ground
------------------------------------------------------------

require "example_base_config"

config = ExampleBaseConfig:new() 

config.fm.atmosphere.ground.creator = ConfigCommon.ground_coxmunk_plus_lamb

config:do_config()

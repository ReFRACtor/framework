------------------------------------------------------------
--- Do a run but with Rayleigh only and no Aerosol.
---
--- Also use LRad, but not LSI
------------------------------------------------------------

require "example_base_config"

config = ExampleBaseConfig:new()

--- Use LRAD only
config.fm.rt.nstream = 4
config.fm.rt.creator = ConfigCommon.radiative_transfer_lrad

--- Use Rayleigh only.
config.fm.atmosphere.aerosol.creator = ConfigCommon.rayleigh_only

config:do_config()

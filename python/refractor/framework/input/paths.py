import os

####

# Paths from environment variables
# Use these to decouple this environment variable from value itself.
# We can change the variables or do additional work here without breaking
# usage of the filenames.

refractor_input_path = os.environ.get('REFRACTOR_INPUT_PATH', './')

#### 

# General data files

aerosol_properties_filename = os.path.join(refractor_input_path, "l2_aerosol_combined.h5")
reference_atmosphere_filename = os.path.join(refractor_input_path, "reference_atmosphere.h5")

####

# Cross section related paths

cross_section_base_path = os.path.join(refractor_input_path, 'cross_sections')

# Maps gas species names to cross section filenames
cross_section_filenames = {
    'O3':     "o3abs_brion_195_660_vacfinal.dat", 
    'H2O':    "H2O_280K900mb.dat",
    'NO2':    "no2r_97.nm",
    'SO2':    "SO2_298_BISA.nm",
    'HCHO':   "H2CO_Meller_Moortgat_MPI.txt",
    'BrO':    "228kbro10cm_padded.nm",
    'CHOCHO': "CHOCHO_Xsections_250-510nm.dat",
}

# Add path to all filenames at once
for gas_name, base_name in cross_section_filenames.items():
    cross_section_filenames[gas_name] = os.path.join(cross_section_base_path, base_name)

# Conversion factors for the files above, ideally we should test rolling these conversions into the files themselves
cross_section_file_conversion_factors = {
    "O3":     1e20,
    "CHOCHO": 5e18,
}

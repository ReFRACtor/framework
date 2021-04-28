import os

####

# Make this import to decouple this environment variable from value itself
refractor_input_path = os.environ['REFRACTOR_INPUT_PATH']

####

# Cross section related paths

cross_section_base_path = os.path.join(refractor_input_path, 'cross_sections')

# Maps gas species names to cross section filenames
cross_section_filenames = {
    'O3':     "o3abs_brion_195_660_vacfinal.dat", 
    'NO2':    "no2r_97.nm",
    'SO2':    "SO2_298_BISA.nm",
    'HCHO':   "H2CO_Meller_Moortgat_MPI.txt",
    'BrO':    "228kbro10cm_padded.nm",
    'CHOCHO': "CHOCHO_Xsections_250-510nm.dat",
}

# Add path to all filenames at once
for gas_name, base_name in cross_section_filenames.items():
    cross_section_filenames[gas_name] = os.path.join(cross_section_base_path, base_name)

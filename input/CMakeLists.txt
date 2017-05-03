# Install files from this directory into a common subdirectory
set(INPUT_INSTALL_PREFIX etc/full_physics)

# Create version.lua file by evaluating cmake variables and install
configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/version.lua.in"
    "${PROJECT_BINARY_DIR}/version.lua"
)
install(FILES "${PROJECT_BINARY_DIR}/version.lua" DESTINATION ${INPUT_INSTALL_PREFIX}/config)

# Common config and inputs
install(DIRECTORY common/config/ DESTINATION ${INPUT_INSTALL_PREFIX}/config
    FILES_MATCHING PATTERN "*.lua")

install(DIRECTORY common/input/ DESTINATION ${INPUT_INSTALL_PREFIX}/input
    FILES_MATCHING PATTERN "*.h5"
        PATTERN "aerosol" EXCLUDE
        PATTERN "solar" EXCLUDE)

# Instrument specific configs and inputs
foreach(sub_dir fts gosat oco)
    install(DIRECTORY ${sub_dir}/config/ DESTINATION ${INPUT_INSTALL_PREFIX}/config
        FILES_MATCHING PATTERN "*.lua" PATTERN "config*.lua" EXCLUDE)
    install(DIRECTORY ${sub_dir}/config/ DESTINATION ${INPUT_INSTALL_PREFIX}/config/${sub_dir}/config
        FILES_MATCHING PATTERN "config*.lua")

    install(DIRECTORY ${sub_dir}/input/ DESTINATION ${INPUT_INSTALL_PREFIX}/input
        FILES_MATCHING PATTERN "*.h5"
            PATTERN "ascii" EXCLUDE
            PATTERN "solar" EXCLUDE
            PATTERN "old" EXCLUDE)
endforeach()

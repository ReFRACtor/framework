# Locate AER OSS Library

# To locate OSS, we expect one of the following is true:
# 1) OSS_DIR will be set to a directory that contains libs in ./lib or ./lib64
#     and includes in ./include
#
# 2) OSS_LIBRARY_DIR will be set to the directory that contains OSS libs and 
#    OSS_INCLUDE_DIR will be set to the directory that contains OSS includes

FIND_PATH(_OSS_INCLUDE_DIR oss_ir_module.mod
  HINTS ${OSS_INCLUDE_DIR} $ENV{OSS_DIR} ${OSS_DIR}
  PATH_SUFFIXES include
)


FIND_LIBRARY(_OSS_LIBRARY
  NAMES oss
  HINTS ${OSS_LIBRARY_DIR} $ENV{OSS_DIR} ${OSS_DIR}
  PATH_SUFFIXES lib64 lib
)

IF(_OSS_LIBRARY AND _OSS_INCLUDE_DIR)
    SET(OSS_LIBRARIES
        "${_OSS_LIBRARY}" CACHE STRING "AER OSS Libraries")
    SET(OSS_INCLUDE_DIRS 
        "${_OSS_INCLUDE_DIR}" CACHE STRING "AER OSS Includes")    
ENDIF(_OSS_LIBRARY AND _OSS_INCLUDE_DIR)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OSS DEFAULT_MSG OSS_LIBRARIES OSS_INCLUDE_DIRS)

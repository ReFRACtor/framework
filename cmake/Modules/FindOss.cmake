# Locate AER OSS Library

# To locate OSS, we expect one of the following is true:
# 1) OSS_DIR will be set to a directory that contains libs in ./lib or ./lib64
#      and includes in ./include
#    GFORTRAN_DIR will be set to a directory that contains libgfortran in ./lib
#      or ./lib64
#
# 2) OSS_LIBRARY_DIR will be set to the directory that contains OSS libs and 
#    OSS_INCLUDE_DIR will be set to the directory that contains OSS includes
#    GFORTRAN_LIBRARY_DIR will be set to the directory that contains libgfortran
FIND_PATH(_OSS_INCLUDE_DIR oss_ir_module.mod
  HINTS ${OSS_INCLUDE_DIR} $ENV{OSS_DIR} ${OSS_DIR}
  PATH_SUFFIXES include include/darwin_x86_64 include/linux_x86_64
)


FIND_LIBRARY(_OSS_LIBRARY
  NAMES oss
  HINTS ${OSS_LIBRARY_DIR} $ENV{OSS_DIR} ${OSS_DIR}
  PATH_SUFFIXES lib64 lib lib/darwin_x86_64 lib/linux_x86_64
)

FIND_LIBRARY(_OSS_FORTRAN_LIBRARY
  NAMES foss
  HINTS ${OSS_LIBRARY_DIR} $ENV{OSS_DIR} ${OSS_DIR}
  PATH_SUFFIXES lib64 lib lib/darwin_x86_64 lib/linux_x86_64
)

FIND_LIBRARY(_OSS_C_LIBRARY
  NAMES coss
  HINTS ${OSS_LIBRARY_DIR} $ENV{OSS_DIR} ${OSS_DIR}
  PATH_SUFFIXES lib64 lib lib/darwin_x86_64 lib/linux_x86_64
)

FIND_LIBRARY(_GFORTRAN_LIBRARY
  NAMES gfortran
  HINTS {GFORTRAN_LIBRARY_DIR} ${GFORTRAN_DIR} $ENV{GFORTRAN_DIR}
  PATH_SUFFIXES lib64 lib lib/darwin_x86_64 lib/linux_x86_64

)

SET(OSS_LIBRARIES "")
IF(_OSS_LIBRARY AND _OSS_FORTRAN_LIBRARY AND _OSS_C_LIBRARY AND _GFORTRAN_LIBRARY AND _OSS_INCLUDE_DIR)
    list(APPEND OSS_LIBRARIES ${_OSS_LIBRARY})
    list(APPEND OSS_LIBRARIES ${_OSS_FORTRAN_LIBRARY})
    list(APPEND OSS_LIBRARIES ${_OSS_C_LIBRARY})
    list(APPEND OSS_LIBRARIES ${_GFORTRAN_LIBRARY})
    SET(OSS_INCLUDE_DIRS 
        "${_OSS_INCLUDE_DIR}" CACHE STRING "AER OSS Includes")    
ENDIF(_OSS_LIBRARY AND _OSS_FORTRAN_LIBRARY AND _OSS_C_LIBRARY AND _GFORTRAN_LIBRARY AND _OSS_INCLUDE_DIR)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OSS DEFAULT_MSG OSS_LIBRARIES OSS_INCLUDE_DIRS)

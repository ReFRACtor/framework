# Files in the current directory to include in the library
set(DIR_SOURCES
    lua_blitz.cc
    lua_callback.cc
    lua_state.cc
    luabind_object.cc
    register_lua.cc
)

set(TEST_SOURCES
    lua_state_test.cc
    lua_test.cc
)


# Full path already included
set(SWIG_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/lua_state.i
    ${CMAKE_CURRENT_SOURCE_DIR}/luabind_object.i
)
list(APPEND SWIG_INTERFACE_FILES ${SWIG_SOURCES})

# Add current source directory to source file names
foreach(source_name ${DIR_SOURCES})
    set(FP_SOURCE_FILES ${FP_SOURCE_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/${source_name})
endforeach(source_name)

foreach(source_name ${TEST_SOURCES})
    set(FP_TEST_FILES ${FP_TEST_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/${source_name})
endforeach(source_name)

# Expose the list of source files to the parent cmake file
set(FP_SOURCE_FILES ${FP_SOURCE_FILES} PARENT_SCOPE)
set(FP_TEST_FILES ${FP_TEST_FILES} PARENT_SCOPE)
set(SWIG_INTERFACE_FILES ${SWIG_INTERFACE_FILES} PARENT_SCOPE)

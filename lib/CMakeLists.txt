
foreach(sub_dir Implementation Lua Python RegisterOutput Support)
    set(LIBRARY_INCLUDE_DIRS ${LIBRARY_INCLUDE_DIRS} "${PROJECT_SOURCE_DIR}/lib/${sub_dir}")
endforeach(sub_dir)

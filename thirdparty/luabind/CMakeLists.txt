cmake_minimum_required(VERSION 3.2)

include(ExternalProject)

project(luabind)

find_package(Lua 5.2 REQUIRED)

# Using updated and fixed luabind version from:
# https://github.com/Oberon00/luabind
# Commit version used is shown in the filename below
set(LUABIND_URL ${CMAKE_CURRENT_SOURCE_DIR}/luabind-Oberon00-60e576e.tar.gz)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include(luabind.cmake)


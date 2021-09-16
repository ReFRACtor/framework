cmake_minimum_required(VERSION 3.2)

include(ExternalProject)

project(luabind)

find_package(Lua 5.2 REQUIRED)

# May want to revisit this, we may want to have the standard c++-17 for
# anaconda here. This should match what framework is built with, we had
# weird errors when luabind was c++-17 but framework was c++-11.
# Global setting: Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Using updated and fixed luabind version from:
# https://github.com/Oberon00/luabind
# Commit version used is shown in the filename below
set(LUABIND_URL ${CMAKE_CURRENT_SOURCE_DIR}/luabind-Oberon00-c0b9359.tar.gz)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include(luabind.cmake)


# Force C++11 or else build will fail if C++17+ is detected due to not
# being able to find timespec_get
set(CMAKE_CXX_STANDARD 11)

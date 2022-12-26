cmake_minimum_required(VERSION 3.24)
project(Map_drift)

set(CMAKE_CXX_STANDARD 17)
set(SOURCE_LIB extra_tools.cpp CubicSlpine.h)

add_library(extra STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем extra

add_executable(Map_drift main.cpp)
target_link_libraries(Map_drift extra)

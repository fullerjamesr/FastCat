cmake_minimum_required(VERSION 3.6)
project(FastCat)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -ffast-math")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto -fopenmp")

set(SOURCE_FILES main.c fileio.c fileio.h coord.h coord.c cpdb.c cpdb.h atominfo.h saxs.c saxs.c saxs.h atominfo.c argparse.c argparse.h)
add_executable(fastcat ${SOURCE_FILES})
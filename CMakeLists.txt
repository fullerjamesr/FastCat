cmake_minimum_required(VERSION 2.8)
project(FastCat)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

if(MINGW)
    add_definitions(-D__USE_MINGW_ANSI_STDIO=1)
endif()

# Optimizer  and warning flags -- supports Intel compiler or falls back to GNU flags
if(CMAKE_C_COMPILER STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${$CMAKE_C_FLAGS} -xHost")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ipo")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -w3")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -ffast-math")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -pedantic")
endif()

# Enable openmp support
FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
endif()

set(SOURCE_FILES main.c fileio.c fileio.h coord.h coord.c cpdb.c cpdb.h atominfo.h saxs.c saxs.c saxs.h atominfo.c argparse.c argparse.h)
add_executable(fastcat ${SOURCE_FILES})
cmake_minimum_required(VERSION 3.10)
project(nucutils C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(H5_UTILS "h5utils")
set(HTSLIB_UTILS "htsutils")

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find MPI
find_package(MPI REQUIRED)

# Find HDF5
find_package(HDF5 REQUIRED)

# Find HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

# Add executable for embed
add_executable(embed embed.c ${H5_UTILS}.c ${HTSLIB_UTILS}.c)
# Add executable for convert
add_executable(convert convert.c)

# Include directories
target_include_directories(embed PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_DIRS}
)

# Link directories
target_link_directories(embed PRIVATE
    ${HTSLIB_LIBRARY_DIRS}
)

# Link libraries
target_link_libraries(embed PRIVATE
    ${HTSLIB_LIBRARIES}
    ${HDF5_LIBRARIES}
    MPI::MPI_C
    OpenMP::OpenMP_C
)

# Compiler options
target_compile_options(embed PRIVATE -O3)
target_compile_options(convert PRIVATE -O3)

if(APPLE)
    # M1 Mac specific flags
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        target_include_directories(embed PRIVATE /opt/homebrew/opt/libomp/include)
        target_link_directories(embed PRIVATE /opt/homebrew/opt/libomp/lib)
    endif()
endif()

# Installation target
install(TARGETS embed convert
    RUNTIME DESTINATION bin
)

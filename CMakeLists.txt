cmake_minimum_required(VERSION 3.25)
project(minimal_robust_positively_invariant_set CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CGAL_HEADER_ONLY ON)

find_package(PkgConfig REQUIRED)

set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:/home/grilo/Git/qhull/install-x86/lib/pkgconfig")
pkg_check_modules(QHULL REQUIRED qhullcpp_d qhull_rd)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)

find_package(PkgConfig REQUIRED)

pkg_check_modules(GMP gmp)
if(NOT ${GMP_FOUND})
        message(FATAL_ERROR "Couldn't find GMP installed")
else()
        message(STATUS "GMP found and installed")
        message(STATUS "${GMP_LINK_LIBRARIES}")
endif()

pkg_check_modules(MPFR mpfr)
if(NOT ${MPFR_FOUND})
        message(FATAL_ERROR "Couldn't find MPFR installed")
else()
        message(STATUS "MPFR found and installed")
        message(STATUS "${MPFR_LINK_LIBRARIES}")
endif()

set(EXE_NAME minimal_robust_positively_invariant_set)

add_executable(${EXE_NAME} main.cpp ../../../../mnt/Data/Git/qhull/src/libqhullcpp/Qhull.cpp)

target_include_directories(${EXE_NAME} PRIVATE ${QHULL_INCLUDE_DIRS} ../gnuplot-iostream)
target_link_directories(${EXE_NAME} PRIVATE ${QHULL_LIBRARY_DIRS})
target_link_libraries(${EXE_NAME}
        PRIVATE
        ${QHULL_LIBRARIES}
        ${ARMADILLO_LIBRARIES}
        Boost::iostreams
        Boost::system
        Boost::filesystem
        Boost::headers
        ${GMP_LINK_LIBRARIES}
        ${MPFR_LINK_LIBRARIES})
target_compile_options(${EXE_NAME} PUBLIC ${QHULL_CFLAGS_OTHER})

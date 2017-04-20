# ==================================
# CMake makefile for CXBUPS
# 18-07-2013 Jelle de Plaa, SRON
# ==================================

cmake_minimum_required (VERSION 2.6)
project (cxbups)


set(CMAKE_C_FLAGS "-g -O2 -msse4.1")


# ==================================
# Find library dependencies 
# ==================================

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_library(M_LIB m)

find_package(pgplot REQUIRED)
find_library(PGPLOT_LIBRARY cpgplot)
find_package(X11 REQUIRED)
find_library(X11_LIBRARY X11)
find_package(GSL REQUIRED)
find_library(GSL_LIBRARY GSL)


# ==================================
# Create liblehmer.a
# ==================================

add_library(lehmer lehmer.c)


# ==================================
# Add executables
# ==================================

add_executable(cxbups cxbups.c)
target_link_libraries(cxbups lehmer ${M_LIB})
add_dependencies(cxbups lehmer)

add_executable(cxbopt cxbopt.c)
target_link_libraries(cxbopt lehmer ${M_LIB} ${PGPLOT_LIBRARY} ${X11_LIBRARY} )
add_dependencies(cxbopt lehmer)

add_executable(cxbrnd cxbrnd.c)
target_link_libraries(cxbrnd lehmer ${M_LIB} ${GSL_LIBRARIES})
add_dependencies(cxbopt lehmer)



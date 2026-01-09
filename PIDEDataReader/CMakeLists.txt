#----------------------------------------------------------------------------
# PIDEDataReader Module
# This module provides functionality to read and parse the GSI PIDE 
# (Particle Irradiation Data Ensemble) database for cell survival experiments.
# 
# Features:
# - Read PIDE Excel/CSV database with LQ parameters
# - Read raw photon and ion survival data (.dat files)
# - Query experiments by cell line, ion type, LET range, etc.
# - Fit LQ parameters to simulation data
# - Compare simulation results with experimental data
#
# No external dependencies - uses standard C++ only
#----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(PIDEDataReader)

#----------------------------------------------------------------------------
# Locate sources and headers for this module
#----------------------------------------------------------------------------
set(PIDE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PIDE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB PIDE_SOURCES ${PIDE_SRC_DIR}/*.cc)
file(GLOB PIDE_HEADERS ${PIDE_INCLUDE_DIR}/*.hh)

#----------------------------------------------------------------------------
# Create the PIDEDataReader library
#----------------------------------------------------------------------------
add_library(PIDEDataReader SHARED ${PIDE_SOURCES} ${PIDE_HEADERS})

target_include_directories(PIDEDataReader PUBLIC ${PIDE_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Build test executable
#----------------------------------------------------------------------------
option(BUILD_PIDE_TEST "Build PIDEDataReader test executable" ON)

if(BUILD_PIDE_TEST)
    add_executable(test_pide_reader 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_pide_reader.cc)
    
    target_link_libraries(test_pide_reader PIDEDataReader)
    
    target_include_directories(test_pide_reader PRIVATE
        ${PIDE_INCLUDE_DIR})
    
    message(STATUS "PIDEDataReader test executable configured")
endif()

#----------------------------------------------------------------------------
# Install the library
#----------------------------------------------------------------------------
install(TARGETS PIDEDataReader DESTINATION lib)
install(FILES ${PIDE_HEADERS} DESTINATION include)

message(STATUS "PIDEDataReader module configured (no external dependencies)")

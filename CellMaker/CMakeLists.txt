#----------------------------------------------------------------------------
# CellMaker Module
# This module provides cell layout initialization functionality
# No external dependencies - uses standard C++ only
#----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(CellMaker)

#----------------------------------------------------------------------------
# Locate sources and headers for this module
#----------------------------------------------------------------------------
set(CELLMAKER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CELLMAKER_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB CELLMAKER_SOURCES ${CELLMAKER_SRC_DIR}/*.cc)
file(GLOB CELLMAKER_HEADERS ${CELLMAKER_INCLUDE_DIR}/*.hh)

#----------------------------------------------------------------------------
# Create the CellMaker library
#----------------------------------------------------------------------------
add_library(CellMaker SHARED ${CELLMAKER_SOURCES} ${CELLMAKER_HEADERS})

target_include_directories(CellMaker PUBLIC ${CELLMAKER_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Build test executables
#----------------------------------------------------------------------------
option(BUILD_CELLMAKER_TEST "Build CellMaker test executables" ON)

if(BUILD_CELLMAKER_TEST)
    # Test: CellLayoutInitializer comprehensive functionality test (standalone, no Geant4)
    add_executable(test_cellLayoutInitializer 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_cellLayoutInitializer.cc)
    
    target_link_libraries(test_cellLayoutInitializer CellMaker)
    
    target_include_directories(test_cellLayoutInitializer PRIVATE
        ${CELLMAKER_INCLUDE_DIR})
    
    # Note: test_cellLayout.cc requires Geant4 (RADCellSimulation) and is built from main CMakeLists.txt
    
    message(STATUS "CellMaker test executable configured (test_cellLayoutInitializer)")
endif()

#----------------------------------------------------------------------------
# Install the library
#----------------------------------------------------------------------------
install(TARGETS CellMaker DESTINATION lib)
install(FILES ${CELLMAKER_HEADERS} DESTINATION include)

message(STATUS "CellMaker module configured (no external dependencies)")

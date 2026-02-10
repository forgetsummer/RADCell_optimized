#----------------------------------------------------------------------------
# CellStateModelCalibration Module
# 
# Provides tools for calibrating cell state model parameters to experimental
# survival curve data using maximum likelihood estimation.
#
# Calibration Methods:
# 1. Direct Lethal (S1->S3): Simple model, direct cell death from damage
# 2. Repair Mediated (S1->S2->S3): Includes repair state with kinetics
#----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(CellStateModelCalibration)

#----------------------------------------------------------------------------
# Setup include and source directories
#----------------------------------------------------------------------------
set(CALIBRATION_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CALIBRATION_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

#----------------------------------------------------------------------------
# Collect source files
#----------------------------------------------------------------------------
file(GLOB CALIBRATION_SOURCES ${CALIBRATION_SRC_DIR}/*.cc)
file(GLOB CALIBRATION_HEADERS ${CALIBRATION_INCLUDE_DIR}/*.hh)

#----------------------------------------------------------------------------
# Create shared library
#----------------------------------------------------------------------------
add_library(CellStateModelCalibration SHARED 
    ${CALIBRATION_SOURCES} 
    ${CALIBRATION_HEADERS}
)

target_include_directories(CellStateModelCalibration PUBLIC
    ${CALIBRATION_INCLUDE_DIR}
)

# C++11 or later required
set_target_properties(CellStateModelCalibration PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)

#----------------------------------------------------------------------------
# Build test executables
#----------------------------------------------------------------------------
option(BUILD_CALIBRATION_TESTS "Build CellStateModelCalibration test executables" ON)

if(BUILD_CALIBRATION_TESTS)
    # Direct Lethal Model Test (S1->S3 only)
    add_executable(test_direct_lethal_calibration 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_direct_lethal_calibration.cc
    )
    
    target_link_libraries(test_direct_lethal_calibration CellStateModelCalibration)
    
    target_include_directories(test_direct_lethal_calibration PRIVATE
        ${CALIBRATION_INCLUDE_DIR}
    )
    
    set_target_properties(test_direct_lethal_calibration PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
    )
    
    message(STATUS "Direct Lethal Calibration test configured")
    
    # Repair Mediated Model Test (S1->S2->S3)
    add_executable(test_repair_mediated_calibration 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_repair_mediated_calibration.cc
    )
    
    target_link_libraries(test_repair_mediated_calibration CellStateModelCalibration)
    
    target_include_directories(test_repair_mediated_calibration PRIVATE
        ${CALIBRATION_INCLUDE_DIR}
    )
    
    set_target_properties(test_repair_mediated_calibration PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
    )

    message(STATUS "Repair Mediated Calibration test configured")

    # Misrepair Model Test (S1->S2->S3 with k_error stochastic misrepair)
    add_executable(test_misrepair_calibration
        ${CMAKE_CURRENT_SOURCE_DIR}/test_misrepair_calibration.cc
    )

    target_link_libraries(test_misrepair_calibration CellStateModelCalibration)

    target_include_directories(test_misrepair_calibration PRIVATE
        ${CALIBRATION_INCLUDE_DIR}
    )

    set_target_properties(test_misrepair_calibration PROPERTIES
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
    )

    message(STATUS "Misrepair Calibration test configured")
endif()

#----------------------------------------------------------------------------
# Install
#----------------------------------------------------------------------------
install(TARGETS CellStateModelCalibration DESTINATION lib)
install(FILES ${CALIBRATION_HEADERS} DESTINATION include)

message(STATUS "CellStateModelCalibration module configured")

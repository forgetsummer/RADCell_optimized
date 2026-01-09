#----------------------------------------------------------------------------
# PhysicalBioTranslator Module
# This module provides cell state and cell cycle modeling functionality
# It translates physical quantities (DNA damage, radiation dose) to 
# biological responses (cell state transitions, cell cycle progression)
# No external dependencies - uses standard C++ only
#----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(PhysicalBioTranslator)

#----------------------------------------------------------------------------
# Locate sources and headers for this module
#----------------------------------------------------------------------------
set(PBT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PBT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB PBT_SOURCES ${PBT_SRC_DIR}/*.cc)
file(GLOB PBT_HEADERS ${PBT_INCLUDE_DIR}/*.hh)

#----------------------------------------------------------------------------
# Create the PhysicalBioTranslator library
#----------------------------------------------------------------------------
add_library(PhysicalBioTranslator SHARED ${PBT_SOURCES} ${PBT_HEADERS})

target_include_directories(PhysicalBioTranslator PUBLIC ${PBT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Build test executables
#----------------------------------------------------------------------------
option(BUILD_PBT_TEST "Build PhysicalBioTranslator test executables" ON)

if(BUILD_PBT_TEST)
    # Test 1: Cell and CellStateModel functionality test
    add_executable(test_cellStateModel 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_cellStateModel.cc)
    
    target_link_libraries(test_cellStateModel PhysicalBioTranslator)
    
    target_include_directories(test_cellStateModel PRIVATE
        ${PBT_INCLUDE_DIR})
    
    # Test 2: Cell Phase Transition test
    add_executable(test_phaseTransition 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_phaseTransition.cc)
    
    target_link_libraries(test_phaseTransition PhysicalBioTranslator CellMaker)
    
    target_include_directories(test_phaseTransition PRIVATE
        ${PBT_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/CellMaker/include)
    
    # Test 3: Cell State Transition test
    add_executable(test_stateTransition 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_stateTransition.cc)
    
    target_link_libraries(test_stateTransition PhysicalBioTranslator CellMaker)
    
    target_include_directories(test_stateTransition PRIVATE
        ${PBT_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/CellMaker/include)
    
    # Test 4: Phase and State vs Dose test
    add_executable(test_phase_state_dose 
        ${CMAKE_CURRENT_SOURCE_DIR}/test_phase_state_dose.cc)
    
    target_link_libraries(test_phase_state_dose PhysicalBioTranslator CellMaker)
    
    target_include_directories(test_phase_state_dose PRIVATE
        ${PBT_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/CellMaker/include)
    
    # Test 5: Phase and State vs Dose test (Parallel version with OpenMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        add_executable(test_phase_state_dose_parallel 
            ${CMAKE_CURRENT_SOURCE_DIR}/test_phase_state_dose_parallel.cc)
        
        target_link_libraries(test_phase_state_dose_parallel 
            PhysicalBioTranslator 
            CellMaker
            OpenMP::OpenMP_CXX)
        
        target_include_directories(test_phase_state_dose_parallel PRIVATE
            ${PBT_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR}/CellMaker/include)
        
        message(STATUS "OpenMP found - building parallel version of test_phase_state_dose")
        
        # Test 6: Full simulation test (adapted from test.cc)
        add_executable(test_full_simulation 
            ${CMAKE_CURRENT_SOURCE_DIR}/test_full_simulation.cc)
        
        target_link_libraries(test_full_simulation 
            PhysicalBioTranslator 
            CellMaker
            OpenMP::OpenMP_CXX)
        
        target_include_directories(test_full_simulation PRIVATE
            ${PBT_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR}/CellMaker/include)
    else()
        message(STATUS "OpenMP not found - skipping parallel version")
    endif()
    
    message(STATUS "PhysicalBioTranslator test executables configured")
endif()

#----------------------------------------------------------------------------
# Install the library
#----------------------------------------------------------------------------
install(TARGETS PhysicalBioTranslator DESTINATION lib)
install(FILES ${PBT_HEADERS} DESTINATION include)

message(STATUS "PhysicalBioTranslator module configured (no external dependencies)")

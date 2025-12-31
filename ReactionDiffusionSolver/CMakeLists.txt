#----------------------------------------------------------------------------
# ReactionDiffusionSolver Module
# This module provides reaction-diffusion simulation functionality
#----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ReactionDiffusionSolver)

#----------------------------------------------------------------------------
# Locate sources and headers for this module
#----------------------------------------------------------------------------
set(RDS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(RDS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB RDS_SOURCES ${RDS_SRC_DIR}/*.cc)
file(GLOB RDS_HEADERS ${RDS_INCLUDE_DIR}/*.hh)

#----------------------------------------------------------------------------
# Create the ReactionDiffusionSolver library
#----------------------------------------------------------------------------
add_library(ReactionDiffusionSolver SHARED ${RDS_SOURCES} ${RDS_HEADERS})

target_include_directories(ReactionDiffusionSolver PUBLIC ${RDS_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Build test executable if CellMaker is available
#----------------------------------------------------------------------------
option(BUILD_RDS_TEST "Build ReactionDiffusionSolver test executable" ON)

if(BUILD_RDS_TEST)
    # Check if CellMaker target exists (when built as part of main project)
    if(TARGET CellMaker)
        # Test 1: ReactionDiffusionSimulation kernel test
        add_executable(test_reactionDiffusion 
            ${CMAKE_CURRENT_SOURCE_DIR}/test_reactionDiffusionKernel.cc)
        
        target_link_libraries(test_reactionDiffusion 
            ReactionDiffusionSolver 
            CellMaker)
        
        target_include_directories(test_reactionDiffusion PRIVATE
            ${RDS_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR}/CellMaker/include)
        
        # Test 2: DiffusionReactionSolver step-by-step test
        add_executable(test_reactionDiffusionSolver1 
            ${CMAKE_CURRENT_SOURCE_DIR}/test_reactionDiffusionSolver1.cc)
        
        target_link_libraries(test_reactionDiffusionSolver1 
            ReactionDiffusionSolver 
            CellMaker)
        
        target_include_directories(test_reactionDiffusionSolver1 PRIVATE
            ${RDS_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR}/CellMaker/include)
        
        # Test 3: Profiling test
        add_executable(test_reactionDiffusionSolver_profile 
            ${CMAKE_CURRENT_SOURCE_DIR}/test_reactionDiffusionSolver_profile.cc)
        
        target_link_libraries(test_reactionDiffusionSolver_profile 
            ReactionDiffusionSolver 
            CellMaker)
        
        target_include_directories(test_reactionDiffusionSolver_profile PRIVATE
            ${RDS_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR}/CellMaker/include)
        
        # Test 4: Minimal I/O test (only writes final step)
        add_executable(test_reactionDiffusionSolver_noIO 
            ${CMAKE_CURRENT_SOURCE_DIR}/test_reactionDiffusionSolver_noIO.cc)
        
        target_link_libraries(test_reactionDiffusionSolver_noIO 
            ReactionDiffusionSolver 
            CellMaker)
        
        target_include_directories(test_reactionDiffusionSolver_noIO PRIVATE
            ${RDS_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR}/CellMaker/include)
        
        message(STATUS "ReactionDiffusionSolver test executables configured")
    else()
        message(STATUS "CellMaker not found - test executables will not be built")
        message(STATUS "To build the tests, include CellMaker module before ReactionDiffusionSolver")
    endif()
endif()

#----------------------------------------------------------------------------
# Install the library
#----------------------------------------------------------------------------
install(TARGETS ReactionDiffusionSolver DESTINATION lib)
install(FILES ${RDS_HEADERS} DESTINATION include)

message(STATUS "ReactionDiffusionSolver module configured")


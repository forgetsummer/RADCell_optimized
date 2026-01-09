#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(RADCellSimulation)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
# Add G4MULTITHREADED definition if Geant4 has multithreading support
if(Geant4_multithreaded_FOUND)
  add_definitions(-DG4MULTITHREADED)
endif()

#----------------------------------------------------------------------------
# Add modular components (these are built first so they can be linked)
#----------------------------------------------------------------------------
add_subdirectory(CellMaker)
add_subdirectory(ReactionDiffusionSolver)
add_subdirectory(PhysicalBioTranslator)
add_subdirectory(PIDEDataReader)

#----------------------------------------------------------------------------
# Locate sources and headers for the main project
# Note: We exclude files that are now in modules
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/CellMaker/include
                    ${PROJECT_SOURCE_DIR}/ReactionDiffusionSolver/include
                    ${PROJECT_SOURCE_DIR}/PhysicalBioTranslator/include
                    ${PROJECT_SOURCE_DIR}/PIDEDataReader/include
                    ${Geant4_INCLUDE_DIR})

# Get all source files except those in modules
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

message("sources" ${sources})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

FIND_PACKAGE(OpenMP REQUIRED) # find openMP to implement parallel computing for loops
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(radcelllib SHARED ${sources}) # creating a library for radcell
target_link_libraries(radcelllib CellMaker ReactionDiffusionSolver PhysicalBioTranslator)

find_package(PythonLibs REQUIRED) #find python libs for calling python function in C++
include_directories(${PYTHON_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Main test executable (full simulation)
#----------------------------------------------------------------------------
add_executable(test test.cc ${sources} ${headers})
target_link_libraries(test ${Geant4_LIBRARIES})
target_link_libraries(test ${PYTHON_LIBRARIES})
target_link_libraries(test CellMaker ReactionDiffusionSolver PhysicalBioTranslator)

#----------------------------------------------------------------------------
# Cell Layout test executable (tests cell layout with Geant4 visualization)
#----------------------------------------------------------------------------
add_executable(test_cellLayout 
    ${PROJECT_SOURCE_DIR}/CellMaker/test_cellLayout.cc 
    ${sources} ${headers})
target_link_libraries(test_cellLayout ${Geant4_LIBRARIES})
target_link_libraries(test_cellLayout ${PYTHON_LIBRARIES})
target_link_libraries(test_cellLayout CellMaker ReactionDiffusionSolver PhysicalBioTranslator)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory
#
file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac 
                    ${PROJECT_SOURCE_DIR}/*.in
                    ${PROJECT_SOURCE_DIR}/*.C 
                    ${PROJECT_SOURCE_DIR}/*.py
                    ${PROJECT_SOURCE_DIR}/README
                    )

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS test DESTINATION bin)
install(TARGETS radcelllib DESTINATION lib)
install(TARGETS CellMaker DESTINATION lib)
install(TARGETS ReactionDiffusionSolver DESTINATION lib)
install(TARGETS PhysicalBioTranslator DESTINATION lib)
install(TARGETS PIDEDataReader DESTINATION lib)

#----------------------------------------------------------------------------
# Python wrapper subdirectory
#----------------------------------------------------------------------------
add_subdirectory(python_wrapper)

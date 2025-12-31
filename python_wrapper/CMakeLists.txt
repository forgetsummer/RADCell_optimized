# I don't really know why these two lines are needed!
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(PYINTERFACE)

#include(${CMAKE_CURRENT_SOURCE_DIR}/FindSWIG.cmake) # this is slightly modified FindSWIG.cmake file -  had to tune it to work with OSX

FIND_PACKAGE(SWIG REQUIRED)

# INCLUDE(${SWIG_USE_FILE})
INCLUDE(UseSWIG)

FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(PythonInterp)

IF(SWIG_FOUND)
  message("FOUND SWIG")
ELSE(SWIG_FOUND)
    message("could not find SWIG in your system")
ENDIF(SWIG_FOUND)

# Use Geant4 libraries from the parent project
SET(LIBS
   radcelllib   
   ${Geant4_LIBRARIES}
   pthread
   ${PYTHON_LIBRARIES}
)


INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)
 
SET_SOURCE_FILES_PROPERTIES(radcell.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(radcell python radcell.i)

SWIG_LINK_LIBRARIES(radcell ${LIBS})

install(TARGETS _radcell DESTINATION bin )

set(python_files_path ${CMAKE_BINARY_DIR}/python_wrapper)

INSTALL(FILES 
    ${python_files_path}/radcell.py
	DESTINATION bin
	)

cmake_minimum_required(VERSION 3.14)
project(ppInteractions CXX)

message(STATUS "Building ppInteractions...")


# Set default build-type to release to enable performance improvements
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# compilation flags
set(CMAKE_CXX_FLAGS "-std=c++11")


# display compiler information
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER} (version ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C++ compilation flags: ${CMAKE_CXX_FLAGS}")


# find CRPropa
include(cmake/FindCRPropa.cmake)

# libraries provided by CRPropa (kiss, eigen3)
include_directories(include ${CRPropa_INSTALL_PREFIX}/libs/kiss)

# assemble includes
include_directories(include 
	${CMAKE_CURRENT_BINARY_DIR}/include 
	${CMAKE_CURRENT_SOURCE_DIR}/../include
	${CRPropa_INCLUDE_DIR}  
	)

# source files for compilation
set(SOURCE_FILES
	${CMAKE_CURRENT_SOURCE_DIR}/src/DecayChargedPion.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/DecayNeutralPion.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/DecayEtaMeson.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/DecayMuon.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleDecay.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/NucleusNucleusInteraction.cc
	)

# decide on whether to use Python or not
option(ENABLE_PYTHON "Create python library via SWIG" On)
message(STATUS "Build Python bindings via SWIG: ${ENABLE_PYTHON}")

# compile with Python support
if(ENABLE_PYTHON)
	find_package(Python COMPONENTS Interpreter Development)
	message(STATUS "Python: Found!")
	message(STATUS "  version ${Python_VERSION}")
	message(STATUS "  executable: ${Python_EXECUTABLE}")
	message(STATUS "  libraries: ${Python_LIBRARIES}")
	message(STATUS "  headers: ${Python_INCLUDE_DIRS}")
	message(STATUS "  site packages: ${Python_SITELIB}")
	if(Python_Development_FOUND)
		message(STATUS "  development libraries: Found!")
	elseif(Python_Development_FOUND)
		message(STATUS "  development libraries: NOT found!")
	endif(Python_Development_FOUND)

	#  use Python_INSTALL_PACKAGE_DIR if provided; otherwise, install in Python_SITELIB 
	set(Python_INSTALL_PACKAGE_DIR "${Python_SITELIB}" CACHE PATH "folder in which the python package is installed")
	message(STATUS "  package install directory: ${Python_INSTALL_PACKAGE_DIR}")

	# find SWIG
	find_package(SWIG 4.0.2 REQUIRED)

	# only support versions of SWIG higher than 4.0
	if(SWIG_VERSION VERSION_LESS 4.0.2)
		message(STATUS "SWIG version too old (${SWIG_VERSION} < 4.0.2). Code will not be compile with Python support.")
		set(SWIG_FOUND False)
	endif()

	if(SWIG_FOUND)
		message(STATUS "SWIG installation directory: ${SWIG_DIR}")
		message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
	else(SWIG_FOUND)
		message(STATUS "SWIG: NOT FOUND (cannot build with Python bindings)")
	endif(SWIG_FOUND)

	# add Python headers
	include_directories(include ${Python_INCLUDE_DIRS})

	# SWIG plugin
	set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/ppInteractions_wrap.cxx PROPERTIES GENERATED True)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ppInteractions_wrap.cxx
		COMMAND swig ${SWIG_MODE_FLAG} -c++ -python
		-I${CRPropa_SWIG_PATH} 
		-I${CRPropa_SWIG_INTERFACE_FILE} 
		-I${CRPropa_INCLUDE_DIR} ${CRPropa_SWIG_DEFINES} 
		-I${CMAKE_SOURCE_DIR}/include 
		-I${Python_INCLUDE_DIRS} 
		-dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/ppInteractions_wrap.cxx 
		-outdir ${CMAKE_CURRENT_BINARY_DIR} 
		${CMAKE_CURRENT_SOURCE_DIR}/python/ppInteractions.i
		)

	# build interface and source in one library
	add_library(ppInteractions MODULE ${CMAKE_CURRENT_BINARY_DIR}/ppInteractions_wrap.cxx ${SOURCE_FILES})
	target_link_libraries(ppInteractions ${CRPropa_LIBRARY} ${Python_LIBRARIES})

elseif(ENABLE_PYTHON) # compile code to work in C++ only, without Python support
	# build interface and source in one library
	add_library(ppInteractions MODULE ${SOURCE_FILES})
	target_link_libraries(ppInteractions ${CRPropa_LIBRARY})

endif(ENABLE_PYTHON)

set_target_properties(ppInteractions PROPERTIES PREFIX "")
set_target_properties(ppInteractions PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(ppInteractions PROPERTIES OUTPUT_NAME "_ppInteractions")


# option(ENABLE_TESTING "Build tests and enable test target" On)
# if(ENABLE_TESTING)
# 	set(GTEST_INCLUDE_PATH ${CRPropa_INSTALL_PREFIX}/../libs/gtest/include)
# 	include_directories(${GTEST_INCLUDE_PATH})
# 	link_directories(${CRPropa_INSTALL_PREFIX}/libs/gtest)
# 	message(STATUS "Google testing suite (from CRPropa): ${CRPropa_INSTALL_PREFIX}/libs/gtest/include")
# 	enable_testing()
# endif(ENABLE_TESTING)


# Install
add_definitions(-DPPInteractions_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS ppInteractions DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/ppInteractions.i DESTINATION swig_interface)

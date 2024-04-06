# search for Python to help find CRPropa (if Python bindings are available)
find_package(Python COMPONENTS Interpreter Development)


# Determine which swig file to use, depending on CRPropa's installation
# For unknown reasons, if built-in is absent this gives segmentation errors in various systems
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" On) 
if(ENABLE_SWIG_BUILTIN)
	set(CRPropa_SWIG_FILE "crpropa-builtin.i")
	set(SWIG_MODE_FLAG "-builtin")
else(ENABLE_SWIG_BUILTIN)
	set(CRPropa_SWIG_FILE "crpropa.i")
	set(SWIG_MODE_FLAG "")
endif(ENABLE_SWIG_BUILTIN)

# try to find CRPropa through Python (if not, provide flags manually)
if(Python_FOUND)
	if(NOT CRPropa_INSTALL_PREFIX)
		execute_process(COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" install_prefix OUTPUT_VARIABLE CRPropa_INSTALL_PREFIX)
	endif(NOT CRPropa_INSTALL_PREFIX)

	if(NOT CRPropa_SWIG_PATH)
		execute_process(COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE CRPropa_SWIG_PATH)
	endif(NOT CRPropa_SWIG_PATH)
endif(Python_FOUND)
list(APPEND CMAKE_PREFIX_PATH ${CRPropa_INSTALL_PREFIX})

if(NOT (CRPropa_SWIG_PATH STREQUAL "") OR NOT CRPropa_SWIG_PATH)
	find_path(CRPropa_SWIG_PATH crpropa.i
		HINTS 
			share/crpropa/swig_interface
			${CRPropa_INSTALL_PREFIX}/share/crpropa/swig_interface	
			${CRPropa_INSTALL_PREFIX}/../python
			$ENV{CRPropa_DIR}/python
			$ENV{CRPropa_DIR}/build/share/crpropa/swig_interface
	)
	# list(APPEND CMAKE_PREFIX_PATH ${CRPropa_SWIG_PATH})
endif()

# CRPropa header
find_path(CRPropa_INCLUDE_DIR CRPropa.h 
	HINTS 
		crpropa 
		include 
		include/crpropa 
		${CRPropa_INSTALL_PREFIX}/../include
		${CRPropa_INSTALL_PREFIX}/include
		$ENV{CRPropa_DIR}/include
		$ENV{CRPropa_DIR}/build
		$ENV{CRPropa_DIR}/build/include
	)

# CRPropa library
find_library(CRPropa_LIBRARY NAMES crpropa libcrpropa
	HINTS 
		crpropa
		lib/crpropa 
		crpropa/lib 
		${CRPropa_INSTALL_PREFIX}/lib
		$ENV{CRPropa_DIR}/build
		$ENV{CRPropa_DIR}/build/lib
	)


# SWIG file for CRPropa
set(CRPropa_SWIG_INTERFACE_FILE "${CRPropa_SWIG_PATH}/${CRPropa_SWIG_FILE}")


if(NOT (${CRPropa_SWIG_INTERFACE_FILE} STREQUAL "") AND NOT (${CRPropa_INCLUDE_DIR} NOT STREQUAL "") AND NOT (${CRPropa_LIBRARY}  STREQUAL ""))
	set(CRPropa_FOUND True)
	message(STATUS "CRPropa install prefix: ${CRPropa_INSTALL_PREFIX}")
	message(STATUS "CRPropa SWIG interface file: ${CRPropa_SWIG_INTERFACE_FILE}")
	message(STATUS "CRPropa include path: ${CRPropa_INCLUDE_DIR}")
	message(STATUS "CRPropa library: ${CRPropa_LIBRARY}")
else()
	set(CRPropa_FOUND False)
	set(CRPropa_HELPER_INSTALL_PREFIX ". CRPropa_INSTALL_PREFIX (path to where CRPropa was installed); usually this flag alone is sufficient to find everything else.")
	set(CRPropa_HELPER_INCLUDE_DIR ". CRPropa_INCLUDE_DIR (directory where CRPropa.h is located; usually CRPropa_INSTALL_PREFIX/include)")
	set(CRPropa_HELPER_LIBRARY ". CRPropa_LIBRARY (path to file libcrpropa.so or libcrpropa.dylib [in OSX])")
	set(CRPropa_HELPER_SWIG_PATH ". CRPropa_SWIG_PATH (path to folder where crpropa.i and crpropa-builtin.i are located; usually usually CRPropa_INSTALL_PREFIX/share/crpropa/swig_interface)")
	message(STATUS 	"CRPropa: NOT FOUND. Try setting the relevant flags manually: \n\t${CRPropa_HELPER_INSTALL_PREFIX} \n\t${CRPropa_HELPER_LIBRARY} \n\t${CRPropa_HELPER_INCLUDE_DIR} \n\t${CRPropa_HELPER_SWIG_PATH}")
	return()
endif()



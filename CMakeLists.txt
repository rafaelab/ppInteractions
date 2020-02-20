cmake_minimum_required(VERSION 2.8.12)

# Change name of the library to something meaningful!
set(PLUGINNAME "NucleusNucleusInteraction")
project(${PLUGINNAME} C CXX)

find_package(PythonInterp)
find_package(PythonLibs)
find_package(SWIG REQUIRED)

# get information on CRPropa
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE CRPROPA_SWIG_PATH)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" install_prefix OUTPUT_VARIABLE CRPROPA_INSTALL_PREFIX)
message(STATUS "CRPropa includes : ${CRPROPA_INSTALL_PREFIX}")

find_path(CRPROPA_SWIG_INTERFACE_FILES crpropa.i HINTS crpropa/ share/crpropa python/ ${CRPROPA_SWIG_PATH})

find_path(CRPROPA_INCLUDES CRPropa.h HINTS crpropa/ include/crpropa ${CRPROPA_INSTALL_PREFIX}/include/)

find_library(CRPROPA_LIB NAMES crpropa libcrpropa HINTS crpropa/ lib/crpropa crpropa/lib ${CRPROPA_INSTALL_PREFIX}/lib/)

message(STATUS "CRPropa SWIG interface files: ${CRPROPA_SWIG_INTERFACE_FILES}")
message(STATUS "CRPropa includes : ${CRPROPA_INCLUDES}")
message(STATUS "CRPropa lib: ${CRPROPA_LIB}")

# SWIG plugin
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx PROPERTIES GENERATED true )
add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  COMMAND swig -c++ -python ${PY3} -I${CRPROPA_SWIG_INTERFACE_FILES} -I${CRPROPA_INCLUDES} -I${CMAKE_SOURCE_DIR}/include -I${PYTHON_INCLUDE_DIRS} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/python/${PLUGINNAME}.i)

include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${CRPROPA_INCLUDES} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})


# build interface and source in one library
add_library(${PLUGINNAME} MODULE 
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/src/NucleusNucleusInteraction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleDecay.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DecayEtaMeson.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DecayChargedPion.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DecayNeutralPion.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/DecayMuon.cc
  )


target_link_libraries(${PLUGINNAME} ${CRPROPA_LIB} ${PYTHON_LIBRARIES})

set_target_properties(${PLUGINNAME} PROPERTIES PREFIX "")
set_target_properties(${PLUGINNAME} PROPERTIES OUTPUT_NAME "_${PLUGINNAME}")

execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(
sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
# Uncomment these lines for installation of the plugin

# install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}.py" DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
# install(TARGETS ${PLUGINNAME} LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")


#option(ENABLE_TESTING "Build tests and enable test target" OFF)
#if(ENABLE_TESTING)
#  include_directories(${CRPROPA_INSTALL_PREFIX}/include/gtest)
#  link_directories(${CRPROPA_INSTALL_PREFIX}/libs/gtest)
#  message(STATUS "Google testin suite (from CRPropa): ${CRPROPA_INSTALL_PREFIX}/libs/gtest/include")
#  if(APPLE)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
#  endif(APPLE)
#
#  enable_testing()
#  add_executable(testMatrix test/testMatrix.cpp)
#  target_link_libraries(testMatrix gtest gtest_main pthread)
#  add_test(testMatrix testMatrix)
#endif(ENABLE_TESTING)
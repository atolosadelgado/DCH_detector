#==========================================================================
#  AIDA Detector description implementation 
#--------------------------------------------------------------------------
# Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
# All rights reserved.
#
# For the licensing terms see $DD4hepINSTALL/LICENSE.
# For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
#
#==========================================================================
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(DCHv2)

IF(NOT TARGET DD4hep::DDCore)
  find_package ( DD4hep REQUIRED )
  include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )
  include ( ${DD4hep_DIR}/cmake/DD4hepBuild.cmake )
  dd4hep_configure_output()
ENDIF()


dd4hep_configure_output()

dd4hep_add_plugin(DCH SOURCES src/DriftChamber_o1_v02.cpp
  USES DD4hep::DDCore DD4hep::DDRec DD4hep::DDCond ROOT::Core ROOT::Geom ROOT::GenVector ROOT::MathCore)
install(TARGETS DCH LIBRARY DESTINATION lib)

# Gaudi Algorithm for DCHdigi
#
# dependencies
# For some reason, the find_package with DD4hep must be the first one
find_package(DD4hep REQUIRED COMPONENTS DDRec DDG4 DDParsers)
# dd4hep_set_compiler_flags() expects DD4hep_SET_RPATH to be set to ON
# otherwise it will not set the rpath when installing
set(DD4HEP_SET_RPATH ON)
dd4hep_set_compiler_flags()
find_package(ROOT COMPONENTS RIO Tree MathCore MathMore)
find_package(EDM4HEP)
find_package(k4FWCore)
find_package(Gaudi)
find_package(GSL)
#-------------------
if(Gaudi_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
  add_subdirectory(DCHdigi)
endif()

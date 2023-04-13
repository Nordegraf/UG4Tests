# Copyright (c) 2023:  G-CSC, Goethe University Frankfurt
# Author: Niklas Conen
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

cmake_minimum_required(VERSION 3.15)
# Name of your plugin and sources.
set(pluginName	UG4Tests)
set(SOURCES		tests.cpp
                unit_tests/vector_tests.cpp
                regression_tests/laplace.cpp)

set(CMAKE_CXX_STANDARD_BACKUP ${CMAKE_CXX_STANDARD})
set(CMAKE_CXX_STANDARD 14)

################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

set(GTEST_LIBS gtest gtest_main gmock gmock_main)

add_executable(ug4tests ${SOURCES})
target_link_libraries(ug4tests PUBLIC ug4 ConvectionDiffusion SuperLU ${GTEST_LIBS})

set(CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD_BACKUP})


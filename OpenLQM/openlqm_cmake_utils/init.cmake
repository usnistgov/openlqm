# Copyright 2025, Noblis, Inc.
# https://fingerprint.nist.gov/openlqm

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at https://www.apache.org/licenses/LICENSE-2.0.

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# OpenLQM was developed by Noblis, Inc. under contract to the National
# Institute of Standards and Technology in 2024-2025, as a port of "LQMetric,"
# which was developed by Noblis, Inc. under contract to the Federal Bureau of
# Investigation's Criminal Justice Information Services Division in 2012-2014.

macro(update_minimum_required MIN_VERSION)
    if(${MIN_VERSION} VERSION_GREATER ${CMAKE_MINIMUM_REQUIRED_VERSION})
        cmake_minimum_required(VERSION ${MIN_VERSION})
    endif()
endmacro()

update_minimum_required(3.30)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Interprocedural optimization (LTO)
option(USE_LTO OFF)

if(USE_LTO)
	set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF)
endif()

if(NOT APPLE)
	if(CMAKE_SIZEOF_VOID_P EQUAL 8)
		set(ARCH "x64")
	else()
		set(ARCH "x86")
	endif()
else()
	set(ARCH "mac")
endif()

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
if(APPLE)
	set(CMAKE_MACOSX_RPATH ON)
	set(CMAKE_INSTALL_RPATH "@executable_path/../lib;@loader_path/../lib")
	set(CMAKE_INSTALL_NAME_DIR "@rpath")

	set(APPLE_RPATH_ARG_1 "-DCMAKE_MACOSX_RPATH=ON")
	set(APPLE_RPATH_ARG_2 "-DCMAKE_INSTALL_RPATH=@executable_path/../lib$<SEMICOLON>@loader_path/../lib")
	set(APPLE_RPATH_ARG_3 "-DCMAKE_INSTALL_NAME_DIR=@rpath")
	set(APPLE_RPATH_ARG_4 "-DCMAKE_SKIP_RPATH=TRUE")
	set(APPLE_RPATH_ARG_5 "-DCMAKE_SKIP_INSTALL_RPATH=TRUE")
	set(APPLE_RPATH_ARG_6 "-DCMAKE_BUILD_WITH_INSTALL_RPATH=FALSE")
elseif(UNIX)
	set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN")

	set(APPLE_RPATH_ARG_1 "-DCMAKE_BUILD_RPATH_USE_ORIGIN=TRUE")
	set(APPLE_RPATH_ARG_2 "-DCMAKE_INSTALL_RPATH=$ORIGIN/../lib:$ORIGIN")
	set(APPLE_RPATH_ARG_4 "-DCMAKE_SKIP_RPATH=TRUE")
	set(APPLE_RPATH_ARG_5 "-DCMAKE_SKIP_INSTALL_RPATH=TRUE")
	set(APPLE_RPATH_ARG_6 "-DCMAKE_BUILD_WITH_INSTALL_RPATH=FALSE")
endif()

set(LIB_DEP_DIR "${ARCH}/$<CONFIG>")

# Add compiler options
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	add_compile_options(
		-Wall
		-Wextra
		-Wshadow
		-Wnon-virtual-dtor
		-pedantic
		-Wpedantic
		-Wold-style-cast
		-Woverloaded-virtual
		-Wsign-conversion
	)
elseif(MSVC)
	add_compile_options(
		/permissive-
		/W4
		/w14640
		/w14624
		/w14928
	)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	add_compile_options(
		-Wnull-dereference
		-Wlogical-op
		-Wmisleading-indentation
	)
endif()

# Set default symbol visibility to hidden for GCC/Clang/AppleClang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
	set(CMAKE_CXX_VISIBILITY_PRESET hidden)
	set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
endif()

# Enable Hot Reload for MSVC compilers if supported.
if(POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

# Minimum OpenCV version defaults to 4.11.0 due to an intermitten compilation error in earlier versions.
set(MINIMUM_OPENCV_VERSION 4.11.0 CACHE STRING "The minimum version of OpenCV required. This should be no lower than 4.5.0 due to a licensing change that took effect in that version.")

# Utility functions:

# Get all properties that cmake supports
if(NOT CMAKE_PROPERTY_LIST)
    execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
    
    # Convert command output into a CMake list
    string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)
endif()
    
function(print_properties)
    message("CMAKE_PROPERTY_LIST = ${CMAKE_PROPERTY_LIST}")
endfunction()
    
function(print_target_properties target)
    if(NOT TARGET ${target})
      message(STATUS "There is no target named '${target}'")
      return()
    endif()

    foreach(property ${CMAKE_PROPERTY_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" property ${property})

        # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
        if(property STREQUAL "LOCATION" OR property MATCHES "^LOCATION_" OR property MATCHES "_LOCATION$")
            continue()
        endif()

        get_property(was_set TARGET ${target} PROPERTY ${property} SET)
        if(was_set)
            get_target_property(value ${target} ${property})
            message("${target} ${property} = ${value}")
        endif()
    endforeach()
endfunction()
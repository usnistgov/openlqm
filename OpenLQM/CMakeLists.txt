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

cmake_minimum_required(VERSION 3.30)
project("OpenLQM" VERSION 1.0.0 LANGUAGES CXX)
include(openlqm_cmake_utils/init.cmake)

# Search file at FILE_PATH for a variable definition like "<VAR_NAME>=xxx" and define that variable in the caller's scope
# if present. Otherwise, continue without error
function(extract_and_set VAR_NAME FILE_PATH)
	execute_process(
		COMMAND ${CMAKE_COMMAND}
			-DVAR_NAME=${VAR_NAME}
			-DFILE_PATH=${FILE_PATH}
			-P "${CMAKE_CURRENT_LIST_DIR}/openlqm_cmake_utils/extract_and_set.cmake"
		ERROR_QUIET
		OUTPUT_VARIABLE EXTRACT_AND_SET_OUTPUT
		OUTPUT_STRIP_TRAILING_WHITESPACE
		RESULT_VARIABLE EXTRACT_AND_SET_RESULT
	)

	if(EXTRACT_AND_SET_OUTPUT MATCHES "${VAR_NAME}=(.*)")
		set(${VAR_NAME} "${CMAKE_MATCH_1}" PARENT_SCOPE)
	else()
		set(${VAR_NAME} "" PARENT_SCOPE)
	endif()
endfunction()

if(STATICALLY_LINK_LIBRARY)
	set(OpenCV_STATIC ON)
endif()

string(REPLACE "," ";" OPENCV_REQUIRED_COMPONENTS "${OPENCV_REQUIRED_COMPONENTS}")
string(REPLACE "," ";" OPENCV_REQUIRED_COMPONENT_TARGETS "${OPENCV_REQUIRED_COMPONENT_TARGETS}")

message(STATUS "VCPKG_INSTALLED_DIR=${VCPKG_INSTALLED_DIR}")
find_package(OpenCV ${MINIMUM_OPENCV_VERSION} REQUIRED COMPONENTS ${OPENCV_REQUIRED_COMPONENTS})
string(TOUPPER ${CMAKE_BUILD_TYPE} build_type_all_caps)
set(IMPORTED_LOCATION_CONFIG_VAR "IMPORTED_LOCATION_${build_type_all_caps}")

# Project-specific convenience variables
set(openlqm_include_dir ${PROJECT_SOURCE_DIR}/include)
configure_file(
	${openlqm_include_dir}/openlqm/openlqm.hpp.in
	${openlqm_include_dir}/openlqm/openlqm.hpp
	@ONLY
)
list(APPEND olqm_lib_source
	openlqm.cpp
	minutia.cpp
	img_util.cpp
	clarity.cpp
)
# Add source to this project's shared library
add_library(olqm SHARED ${olqm_lib_source})
set_target_properties(olqm PROPERTIES OUTPUT_NAME openlqm)

# Generate debug info and profiler flags in the Release build.
if(PROFILE_RELEASE_BUILD)
	if(MSVC)
		add_compile_options("$<$<CONFIG:RELEASE>:/Zi>")
		add_link_options("$<$<CONFIG:RELEASE>:/DEBUG>")
		add_link_options("$<$<CONFIG:RELEASE>:/PROFILE>")
	else()
		add_link_options("$<$<CONFIG:RELEASE>:-g>")
		add_link_options("$<$<CONFIG:RELEASE>:-pg>")
	endif()
endif()

set(OPENLQM_MODEL_PATH ${PROJECT_SOURCE_DIR}/models/openlqm_model.xml)
if(WIN32)
	# For local builds, olqm needs a way to find openlqm_model.xml
	set(OPENLQM_MODEL_INSTALL_PATH ${CMAKE_INSTALL_PREFIX}/bin/openlqm_model.xml)
	set(OPENLQM_MODEL_INSTALL_PATH "\"${OPENLQM_MODEL_INSTALL_PATH}\"")
	target_compile_definitions(olqm PRIVATE OPENLQM_MODEL_INSTALL_PATH=${OPENLQM_MODEL_INSTALL_PATH})
endif()

target_include_directories(olqm PUBLIC
    $<BUILD_INTERFACE:${openlqm_include_dir}>
    $<INSTALL_INTERFACE:include>
)

# In order for a relocatable installation of olqm to be able to find itself when loading
# its openlqm_model.xml, it needs to be compiled with POSITION_INDEPENDENT_CODE in order for dladdr
# to be guaranteed to return the shared library's location instead of the location of
# the executable linking to it (which will NOT be OpenLQM package's location when a
# third-arty executable is linking to olqm)
if(NOT WIN32)
	set_target_properties(olqm PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if(STATICALLY_LINK_LIBRARY)
	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
		target_link_options(olqm PRIVATE -static-libstdc++ -static-libgcc)
	elseif(MSVC)
		set_target_properties(olqm PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

target_compile_definitions(olqm PRIVATE OPENLQM_EXPORTS)

# Link opencv to the OpenLQM libary
target_link_libraries(olqm PRIVATE ${OpenCV_LIBS})

# Add source to this project's executable.
add_executable(OpenLQM main.cpp)
set_target_properties(OpenLQM PROPERTIES OUTPUT_NAME openlqm)
target_include_directories(OpenLQM PRIVATE "${openlqm_include_dir}")

# Link the executable to the OpenLQM library
target_link_libraries(OpenLQM PRIVATE olqm)

if(STATICALLY_LINK_EXECUTABLE)
	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
		target_link_options(OpenLQM PRIVATE -static-libstdc++ -static-libgcc) # Clang supports -static-libgcc for compatibility, but -static-libclang_rt can be used more specifically
	elseif(MSVC)
		set_target_properties(OpenLQM PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

# On Windows and Mac, don't split packages
set(runtime_component all)
set(dev_component all)
if(NOT APPLE)
	if(UNIX)
		set(runtime_component Runtime)
		set(dev_component Development)
	endif()
endif()

install(TARGETS OpenLQM
	EXPORT OpenLQMTargets
	COMPONENT ${runtime_component}
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

install(TARGETS olqm
	EXPORT OpenLQMTargets
	COMPONENT ${runtime_component}
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)


if(WIN32)
	install(
		FILES ${install_target_list}
		DESTINATION bin
		COMPONENT ${runtime_component}
	)
	install(
		FILES ${OPENLQM_MODEL_PATH}
		DESTINATION bin
		COMPONENT ${runtime_component}
	)

	# On Windows, install all needed DLLs
	install(CODE [[
	    file(GET_RUNTIME_DEPENDENCIES
		EXECUTABLES $<TARGET_FILE:OpenLQM>
		RESOLVED_DEPENDENCIES_VAR RDEPS
		UNRESOLVED_DEPENDENCIES_VAR UDEPS
		PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*"
		POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
	    )
	    foreach(DEP IN LISTS RDEPS)
		file(INSTALL "${DEP}" DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
		     FOLLOW_SYMLINK_CHAIN)
	    endforeach()

	    foreach(DEP IN LISTS UDEPS)
	        message(WARNING "Unresolved dependency: ${DEP}")
	    endforeach()
]] COMPONENT ${runtime_component})

else()
	install(
		FILES ${install_target_list}
		DESTINATION lib
		COMPONENT ${runtime_component}
	)
	foreach(opencv_target_name ${opencv_target_names})
		set(real_target_location ${CMAKE_INSTALL_PREFIX}/lib/${opencv_target_name})
		if(APPLE)
			install(CODE "
				message(STATUS \"Patching RPATH for ${real_target_location}\")
				execute_process(
						COMMAND install_name_tool -add_rpath \"@executable_path/../lib\" ${real_target_location}
				)
				execute_process(
						COMMAND install_name_tool -add_rpath \"@loader_path/../lib\" ${real_target_location}
				)
			" COMPONENT ${runtime_component})
		else()
			install(CODE "
				message(STATUS \"Removing RUNPATH for ${real_target_location}\")
				execute_process(
						COMMAND chrpath -d \"${real_target_location}\"
				)
			" COMPONENT ${runtime_component})
		endif()
	endforeach()
	install(
		FILES ${OPENLQM_MODEL_PATH}
		DESTINATION share
		COMPONENT ${runtime_component}
	)
endif()
# API header
install(
	FILES include/openlqm/openlqm.hpp
	DESTINATION include/openlqm
	COMPONENT ${dev_component}
)

include(CMakePackageConfigHelpers)
set(OPENLQM_CONFIG_REL_INSTALL_PATH lib/cmake/OpenLQM)
install(EXPORT OpenLQMTargets DESTINATION ${OPENLQM_CONFIG_REL_INSTALL_PATH} COMPONENT ${dev_component})
configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/OpenLQMConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/OpenLQMConfig.cmake
	INSTALL_DESTINATION ${OPENLQM_CONFIG_REL_INSTALL_PATH}
)
write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/OpenLQMConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY AnyNewerVersion
)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/OpenLQMConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/OpenLQMConfigVersion.cmake
	DESTINATION ${OPENLQM_CONFIG_REL_INSTALL_PATH}
	COMPONENT ${dev_component}
)

if(NOT WIN32)
	install(FILES
		${CMAKE_CURRENT_SOURCE_DIR}/openlqm.1
		DESTINATION share/man/man1
		COMPONENT ${runtime_component}
	)
endif()

set(CPACK_PACKAGE_NAME OpenLQM)
set(CPACK_PACKAGE_VENDOR NIST)
set(CPACK_PACKAGE_VERSION 1.0.0)
set(CPACK_PACKAGE_CONTACT OpenLQM@list.nist.gov)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenLQM - A C++ library and CLI for fingerprint analysis")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://fingerprint.nist.gov/openlqm")
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.txt)

if(WIN32)
	find_program(WIX_EXECUTABLE wix)
	set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenLQM")
	if(WIX_EXECUTABLE)
		set(CPACK_WIX_VERSION 4)
		list(APPEND CPACK_GENERATOR "WIX")
		list(APPEND CPACK_SOURCE_GENERATOR "WIX")

		set(CPACK_WIX_UPGRADE_GUID "4275162D-A380-4757-8310-4C2409270811")
		set(CPACK_WIX_PATCH_FILE ${CMAKE_CURRENT_SOURCE_DIR}/openlqm_cmake_utils/package.wxs)
	else()
		set(CPACK_GENERATOR NSIS)
		set(CPACK_SOURCE_GENERATOR NSIS)
		set(CPACK_NSIS_DISPLAY_NAME "${CMAKE_PACKAGE_NAME} ${CMAKE_PACKAGE_VERSION}")
		set(CPACK_NSIS_URL_INFO_ABOUT ${CPACK_PACKAGE_HOMEPAGE_URL})
		set(CPACK_NSIS_CONTACT ${CPACK_PACKAGE_CONTACT})
		set(CPACK_NSIS_MODIFY_PATH ON)
	endif()
elseif(APPLE)
	set(CPACK_GENERATOR "DragNDrop;TGZ;productbuild")
	set(CPACK_SOURCE_GENERATOR "DragNDrop;TGZ;productbuild")
	set(CPACK_COMPONENTS_ALL all)
	set(CPACK_COMPONENT_ALL_DISPLAY_NAME All)
	set(CPACK_PRODUCTBUILD_DOMAINS TRUE)
	set(CPACK_PRODUCTBUILD_DOMAINS_ROOT TRUE)

	execute_process(COMMAND sw_vers -productVersion OUTPUT_VARIABLE OSX_VERS)
	string(STRIP "${OSX_VERS}" OSX_VERS)
	execute_process(COMMAND uname -m OUTPUT_VARIABLE MAC_ARCH)
	string(STRIP "${MAC_ARCH}" MAC_ARCH)
	set(CPACK_PACKAGE_FILE_NAME "openlqm-${CPACK_PACKAGE_VERSION}-macos-${OSX_VERS}-${MAC_ARCH}")
	set(CPACK_PRODUCTBUILD_IDENTIFIER "gov.nist.openlqm")

	# Must be manually set to override ~/Applications and/or /Users/Shared/Applications
	set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local/openlqm)
else()
	set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local/openlqm)
	list(APPEND LINUX_GENERATORS "TGZ")
	find_program(DPKG_DEB_EXECUTABLE dpkg-deb)
	if(DPKG_DEB_EXECUTABLE)
		list(APPEND LINUX_GENERATORS "DEB")
		# We're customizing the package name, so we need to define these ourselves
		set(CPACK_DEBIAN_PACKAGE_RELEASE 1)
		execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
		string(STRIP ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE} CPACK_DEBIAN_PACKAGE_ARCHITECTURE)

		set(CPACK_DEBIAN_DEVELOPMENT_FILE_NAME "openlqm-dev_${CPACK_PACKAGE_VERSION}-${CPACK_DEBIAN_PACKAGE_RELEASE}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}.deb")
		set(CPACK_DEBIAN_RUNTIME_FILE_NAME "openlqm_${CPACK_PACKAGE_VERSION}-${CPACK_DEBIAN_PACKAGE_RELEASE}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}.deb")
		set(CPACK_DEBIAN_DEVELOPMENT_PACKAGE_NAME "openlqm-dev")
		set(CPACK_DEBIAN_RUNTIME_PACKAGE_NAME "openlqm")
	endif()

    find_program(RPMBUILD_EXECUTABLE rpmbuild)
    if(RPMBUILD_EXECUTABLE)
		list(APPEND LINUX_GENERATORS "RPM")
	endif()

	set(CPACK_GENERATOR "${LINUX_GENERATORS}")
	set(CPACK_SOURCE_GENERATOR "${LINUX_GENERATORS}")

	set(CPACK_DEBIAN_PACKAGE_MAINTAINER ${CPACK_PACKAGE_CONTACT})
	set(CPACK_DEB_COMPONENT_INSTALL ON)
	set(CPACK_RPM_PACKAGE_LICENSE "Apache-2.0")
endif()

# Move repo documentation into usable form for packages
find_program(PANDOC pandoc)
if (PANDOC)
	set(DOC_FORMAT "rtf")

	# Strip excess things from README
	execute_process(COMMAND grep -v ^\\[! "${CMAKE_CURRENT_SOURCE_DIR}/../README.md" OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/README.md")
	execute_process(COMMAND sed -i.bak -e s|<img.*/>|| "${CMAKE_CURRENT_BINARY_DIR}/README.md")

	execute_process(COMMAND ${PANDOC} -f markdown_mmd -s --columns 1 -o "${CMAKE_CURRENT_BINARY_DIR}/README.${DOC_FORMAT}" "${CMAKE_CURRENT_BINARY_DIR}/README.md")
	set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_BINARY_DIR}/README.${DOC_FORMAT}")
else()
	file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../README.md" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
	file(RENAME "${CMAKE_CURRENT_BINARY_DIR}/README.md" "${CMAKE_CURRENT_BINARY_DIR}/README.txt")
	set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_BINARY_DIR}/README.txt")
endif()

include(CPack)

if(NOT APPLE)
	if(UNIX)
		cpack_add_component(
			DESCRIPTION "The executable, library, dependencies, and model"
		)
		cpack_add_component(
			DESCRIPTION "Headers and CMake files"
			DEPENDS Runtime
		)
	endif()
endif()

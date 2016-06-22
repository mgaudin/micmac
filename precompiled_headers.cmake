#
# Christoph Heindl 2010
# Precompiled Headers Demo
# http://cheind.wordpress.com
#

# Instructs the MSVC toolset to use the precompiled header PRECOMPILED_HEADER
# for each source file given in the collection named by SOURCE_VARIABLE_NAME.
if(WITH_HEADER_PRECOMP)
	if (MSVC)
		function(enable_precompiled_headers_msvc PRECOMPILED_HEADER SOURCE_VARIABLE_NAME)
			set(files ${${SOURCE_VARIABLE_NAME}})

			# Generate precompiled header translation unit
			get_filename_component(pch_basename ${PRECOMPILED_HEADER} NAME_WE)
			set(pch_abs ${PRECOMPILED_HEADER})
			set(pch_unity ${pch_basename}.cpp)
			FILE(WRITE ${pch_unity} "// Precompiled header unity generated by CMake\n")
			FILE(APPEND ${pch_unity} "#include <${pch_abs}>\n")
			set_source_files_properties(${pch_unity}  PROPERTIES COMPILE_FLAGS "/Yc\"${pch_abs}\"")

			# Update properties of source files to use the precompiled header.
			# Additionally, force the inclusion of the precompiled header at beginning of each source file.
			foreach(source_file ${files} )
			set_source_files_properties(
			${source_file} 
			PROPERTIES COMPILE_FLAGS
			"/Yu\"${pch_abs}\" /FI\"${pch_abs}\""
			)
			endforeach(source_file)

			# Finally, update the source file collection to contain the precompiled header translation unit
			set(${SOURCE_VARIABLE_NAME} ${${SOURCE_VARIABLE_NAME}} ${pch_unity} PARENT_SCOPE)
		endfunction(enable_precompiled_headers_msvc)
	elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_IS_CLANG)
		function(enable_precompiled_headers_GCC HEADER_NAME TARGET_NAME EXTRA_CXX_FLAGS)
			GET_FILENAME_COMPONENT(_name ${HEADER_NAME} NAME)
			SET(_source "${PROJECT_SOURCE_DIR}/include/${HEADER_NAME}")

			if(CMAKE_CXX_COMPILER_IS_CLANG)
				set(precompiled_header "${CMAKE_SOURCE_DIR}/include/${HEADER_NAME}.pch")
				set(pch_flags -x c++-header ${_source} -o ${precompiled_header})
				set(pch_include_flag "-include-pch ${precompiled_header}")
			elseif(CMAKE_COMPILER_IS_GNUCXX)
				set(precompiled_header "${CMAKE_SOURCE_DIR}/include/${HEADER_NAME}.gch")
				set(pch_flags ${_source} -o ${precompiled_header})
				set(pch_include_flag "-include ${HEADER_NAME}")
			endif()
			set(precompiled_header ${precompiled_header} PARENT_SCOPE)
			SET(OPTION_HP "-x c++-header")

			STRING(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" _flags_var_name)
			SET(_compiler_FLAGS ${${_flags_var_name}} ${EXTRA_CXX_FLAGS} ${CMAKE_CXX_FLAGS})

			GET_DIRECTORY_PROPERTY(_directory_flags INCLUDE_DIRECTORIES)
			include_dirs_to_flags("${_directory_flags}" directory_includes)
			list(REMOVE_DUPLICATES directory_includes)
			list(APPEND _compiler_FLAGS ${directory_includes})

			GET_DIRECTORY_PROPERTY(_directory_flags COMPILE_DEFINITIONS)
			foreach(item ${_directory_flags})
				LIST(APPEND _compiler_FLAGS "-D${item}")
			endforeach()

			GET_DIRECTORY_PROPERTY(_directory_flags COMPILE_OPTIONS)
			LIST(APPEND _compiler_FLAGS ${_directory_flags})

			list(REMOVE_DUPLICATES _compiler_FLAGS)
			SEPARATE_ARGUMENTS(_compiler_FLAGS)

			ADD_CUSTOM_COMMAND(
				OUTPUT ${precompiled_header}
				COMMAND ${CMAKE_CXX_COMPILER} ${_compiler_FLAGS} ${pch_flags} -Winvalid-pch
				DEPENDS ${_source} IMPLICIT_DEPENDS CXX ${_source})
				ADD_CUSTOM_TARGET(${TARGET_NAME}_${EXT_HP} DEPENDS ${precompiled_header})

			ADD_DEPENDENCIES(${TARGET_NAME} ${TARGET_NAME}_${EXT_HP})

			foreach(file ${Elise_Src_Files})
				if (NOT ${file} MATCHES "tiff/el_dcraw.c" AND
					NOT ${file} MATCHES "mullgesuhlig/muvmblock.cpp" AND
					NOT ${file} MATCHES "mullgesuhlig/mubasic.cpp" AND
					NOT ${file} MATCHES "mullgesuhlig/muflaguer.cpp" AND
					NOT ${file} MATCHES "mullgesuhlig/mufmueller.cpp")
					set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "${pch_include_flag}")
				endif()
			endforeach()
		endfunction(enable_precompiled_headers_GCC)
	endif()
endif()

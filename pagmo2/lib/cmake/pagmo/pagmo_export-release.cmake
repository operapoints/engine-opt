#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Pagmo::pagmo" for configuration "Release"
set_property(TARGET Pagmo::pagmo APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Pagmo::pagmo PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "TBB::tbb"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libpagmo.so.8.0"
  IMPORTED_SONAME_RELEASE "libpagmo.so.8"
  )

list(APPEND _cmake_import_check_targets Pagmo::pagmo )
list(APPEND _cmake_import_check_files_for_Pagmo::pagmo "${_IMPORT_PREFIX}/lib/libpagmo.so.8.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

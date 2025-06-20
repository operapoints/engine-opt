# Get current dir.
get_filename_component(_PAGMO_CONFIG_SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Find the deps. Alter the cmake module path.
set(_PAGMO_CONFIG_OLD_MODULE_PATH "${CMAKE_MODULE_PATH}")
list(APPEND CMAKE_MODULE_PATH "${_PAGMO_CONFIG_SELF_DIR}")
set(THREADS_PREFER_PTHREAD_FLAG YES)
find_package(Threads REQUIRED)
unset(THREADS_PREFER_PTHREAD_FLAG)
include(PagmoFindBoost)

# Restore original module path.
set(CMAKE_MODULE_PATH "${_PAGMO_CONFIG_OLD_MODULE_PATH}")
unset(_PAGMO_CONFIG_OLD_MODULE_PATH)

include(${_PAGMO_CONFIG_SELF_DIR}/pagmo_export.cmake)

# Clean up.
unset(_PAGMO_CONFIG_SELF_DIR)

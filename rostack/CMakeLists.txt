cmake_minimum_required(VERSION 3.15)

find_package(Pagmo REQUIRED)

add_library(rostack
    src/rostack_utils.cpp
    src/rostack_physics.cpp
)

target_include_directories(rostack
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Also applies C++20 to executables linking against rostack.
target_compile_features(rostack PUBLIC cxx_std_20)

add_executable(rostack_corotating
    src/demos/rostack_corotating.cpp
)
target_link_libraries(rostack_corotating
    PRIVATE rostack Pagmo::pagmo
)

add_executable(rostack_counterrotating
    src/demos/rostack_counterrotating.cpp
)
target_link_libraries(rostack_counterrotating
    PRIVATE rostack Pagmo::pagmo
)
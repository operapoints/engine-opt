cmake_minimum_required(VERSION 3.15)

find_package(Pagmo REQUIRED)

add_library(rostack
    src/rostack_utils.cpp
    src/study_physics.cpp
    src/rostack_optimizer.cpp
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

add_executable(optimize_stack
    src/demos/optimize_stack.cpp
)
target_link_libraries(optimize_stack
    PRIVATE rostack Pagmo::pagmo
)

add_executable(reproduce_python_optimizer
    src/demos/reproduce_python_optimizer.cpp
)
target_link_libraries(reproduce_python_optimizer
    PRIVATE rostack Pagmo::pagmo
)
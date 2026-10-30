cmake_minimum_required(VERSION 3.15)



add_library(rostack
    src/rostack_utils.cpp
    src/rostack_physics.cpp
)

target_include_directories(rostack PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


add_executable(rostack_test src/demos/rostack_test.cpp)
target_link_libraries(rostack_test PRIVATE rostack)

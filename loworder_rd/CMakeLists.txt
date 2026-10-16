cmake_minimum_required(VERSION 3.15)

find_package(Eigen3 REQUIRED NO_MODULE)


add_library(LORd 
src/hessian_H0.cpp
src/solvers.cpp
src/LORd.cpp)

target_include_directories(LORd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(LORd PUBLIC Eigen3::Eigen)
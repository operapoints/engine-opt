cmake_minimum_required(VERSION 3.15)

find_package(Eigen3 REQUIRED NO_MODULE)


add_library(LORd 
src/hessian_H0.cpp
src/solvers.cpp
src/LORd.cpp
src/ASTRA_hessian_H0.cpp)

target_include_directories(LORd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(LORd PUBLIC Eigen3::Eigen)


add_executable(LORd_test src/LORd_demos/LORd_test.cpp)
target_link_libraries(LORd_test PRIVATE LORd)


add_executable(campbell src/LORd_demos/campbell.cpp)   # plus your LORd sources if they aren't a library
target_include_directories(campbell PRIVATE
    ${CMAKE_SOURCE_DIR}/external/matplotlib-cpp   # directory containing matplotlibcpp.h
)
target_link_libraries(campbell PRIVATE
    Python3::Python
    Python3::NumPy
    LORd
)


#add_executable(ASTRA_Hessian_test src/ASTRA_hessian_H0_test.cpp)
#target_link_libraries(ASTRA_Hessian_test PUBLIC Eigen3::Eigen)
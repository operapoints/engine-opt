cmake_minimum_required(VERSION 3.27)

project(Engine_Opt LANGUAGES CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG " -g -O0")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O3 -DNDEBUG")
#list(APPEND CMAKE_PREFIX_PATH "./pagmo2")
find_package(Pagmo REQUIRED)
if(Pagmo_FOUND)
    message(STATUS "Found pagmo at ${Pagmo_DIR}")
endif()

find_package(Boost REQUIRED)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
message(STATUS "Python3 exe:      ${Python3_EXECUTABLE}")
message(STATUS "Python3 library:  ${Python3_LIBRARIES}")
message(STATUS "NumPy include:    ${Python3_NUMPY_INCLUDE_DIRS}")

add_executable(engine_opt
    engine_opt.cpp
    jet_calc.cpp
    jet_calc.h)


target_link_libraries(engine_opt PRIVATE Pagmo::pagmo)

target_include_directories(engine_opt PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_options(engine_opt PRIVATE -lineinfo)

add_subdirectory(loworder_rd)
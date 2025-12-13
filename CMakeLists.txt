cmake_minimum_required(VERSION 3.27)

project(Engine_Opt LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG " -g -O0")
# string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O3 -DNDEBUG")
list(APPEND CMAKE_PREFIX_PATH "./pagmo2")
find_package(Pagmo REQUIRED)

find_package(Boost REQUIRED)


add_executable(engine_opt
    engine_opt.cpp
    jet_calc.cpp
    jet_calc.h)


target_link_libraries(engine_opt PRIVATE Pagmo::pagmo)

target_include_directories(engine_opt PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_options(engine_opt PRIVATE -lineinfo)


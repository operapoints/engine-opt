cmake_minimum_required(VERSION 3.27)

# set(CMAKE_PREFIX_PATH "~/spack/opt/spack/linux-ubuntu24.04-skylake/gcc-13.3.0/pagmo2-2.18.0-xprq2kk2yif7hbppl66uawyg45h6qrhb")

# The name of our project.
project(Engine_Opt LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Debug)

# Look for an installation of pagmo in the system.
find_package(Pagmo REQUIRED)

find_package(Boost REQUIRED)

# Create an executable, and link it to the Pagmo::pagmo imported target.
# This ensures that, in the compilation of 'getting_started', pagmo's include
# dirs are added to the include path of the compiler and that pagmo's
# dependencies are transitively linked to 'getting_started'.
add_executable(engine_opt
    engine_opt.cpp
    jet_calc.cpp
    jet_calc.h)
target_link_libraries(engine_opt PRIVATE Pagmo::pagmo)

target_include_directories(engine_opt PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_options(engine_opt PRIVATE -lineinfo)
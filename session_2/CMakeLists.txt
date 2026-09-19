cmake_minimum_required(VERSION 3.16)
project(kinetic_fisher LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(HDF5 REQUIRED COMPONENTS CXX)

add_library(pic_core STATIC
    src/boris_pusher.cpp
    src/step1_single_particle.cpp
    src/interpolation.cpp
    src/step2_pic_loop.cpp
)
target_compile_options(pic_core PRIVATE -Wall -Wextra)
target_include_directories(pic_core PUBLIC src)

add_executable(kinetic_fisher src/main.cpp)

target_compile_options(kinetic_fisher PRIVATE -Wall -Wextra)
target_include_directories(kinetic_fisher PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(kinetic_fisher PRIVATE pic_core ${HDF5_CXX_LIBRARIES})

enable_testing()

add_executable(test_boris_pusher tests/test_boris_pusher.cpp)
target_link_libraries(test_boris_pusher PRIVATE pic_core)
add_test(NAME test_boris_pusher COMMAND test_boris_pusher)

add_executable(test_deposit_gather tests/test_deposit_gather.cpp)
target_link_libraries(test_deposit_gather PRIVATE pic_core)
add_test(NAME test_deposit_gather COMMAND test_deposit_gather)

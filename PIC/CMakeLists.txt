cmake_minimum_required(VERSION 3.16)
project(pic LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release = -O3 optimisation.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# The physics, plain C++ with no Python in it. Position-independent code because
# it is linked into the Python module, which is a shared library.
add_library(pic_core STATIC
    src/field_solver.cpp
    src/grid.cpp
    src/hybrid.cpp
    src/initial_conditions.cpp
    src/population.cpp
)
target_include_directories(pic_core PUBLIC src)
target_compile_options(pic_core PRIVATE -Wall -Wextra)
# OpenMP: the #pragma omp lines of push and deposit share the particles between the cores.
# The number of threads is set at run time: OMP_NUM_THREADS=4 ./hybrid timing
find_package(OpenMP REQUIRED)
target_link_libraries(pic_core PUBLIC OpenMP::OpenMP_CXX)
set_target_properties(pic_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# cmake -S . -B build -DPIC_NATIVE=ON: use every vector instruction of this CPU
# (session 4: -march=native). Off by default so the build works on any machine.
option(PIC_NATIVE "Optimise for the CPU of this machine" OFF)
if(PIC_NATIVE)
    target_compile_options(pic_core PUBLIC -march=native)
endif()

# The program: ./hybrid <case> runs one of the initial conditions and writes <case>.csv.
add_executable(hybrid src/main.cpp)
target_link_libraries(hybrid PRIVATE pic_core)
target_compile_options(hybrid PRIVATE -Wall -Wextra)

# The Python module `pic`, built into build/ (tests/test_pic.py looks for it there).
pybind11_add_module(pic python/bindings.cpp)
target_link_libraries(pic PRIVATE pic_core)

# All tests are in Python: `ctest` just runs them.
enable_testing()
add_test(NAME python_tests
         COMMAND ${Python_EXECUTABLE} -m unittest discover -s ${CMAKE_SOURCE_DIR}/tests -v)

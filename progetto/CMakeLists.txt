cmake_minimum_required(VERSION 3.5)
set(MY_PROJECT_VERSION 1.0.0)
project(pyRA VERSION ${MY_PROJECT_VERSION}
    DESCRIPTION "Python Evolutionary Reaction Automata"
    LANGUAGES CXX)


set(LIB_NAME "reaction_automata")
set(PY_LIB_NAME "pyRA")
set(HEADER "ReactionAutomata.hpp")
set(SRC "ReactionAutomata.cpp")
set(PY_SRC "ReactionAutomata_py.cpp")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include/")
set(MAIN "main.cpp")
set(PATH_TO_PYBIND11 "" CACHE STRING "Specify the path to the pybind11 module location")
list(APPEND CMAKE_PREFIX_PATH "${PATH_TO_PYBIND11}")
set(PYLIB_INSTALL_PATH "${CMAKE_SOURCE_DIR}/my_program/python_module/" CACHE STRING "Specify the path to install the python module")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/my_program/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/my_program/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/my_program/bin")


add_compile_options(-fopenmp -O3)
find_package(pybind11 REQUIRED)

add_subdirectory(src)

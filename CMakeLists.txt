cmake_minimum_required(VERSION 3.15)
project(CMakeProject1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Znajdź pybind11 i matplot++ globalnie
include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.11.1
)

FetchContent_Declare(
  matplotplusplus
  GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus.git
  GIT_TAG        master
)

FetchContent_MakeAvailable(pybind11 matplotplusplus)


# Dodaj podprojekt src/
add_subdirectory(src)

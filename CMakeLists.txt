cmake_minimum_required(VERSION 3.16)
project(vol_smile_pricer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(pricer_lib
    src/black_scholes.cpp
    src/binomial_tree.cpp
)
target_include_directories(pricer_lib PUBLIC include)

enable_testing()
add_subdirectory(tests)
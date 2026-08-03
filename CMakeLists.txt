cmake_minimum_required(VERSION 3.16)
project(vol_smile_pricer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(pricer_lib
    src/black_scholes.cpp
    src/binomial_tree.cpp
    src/delta_hedge.cpp
)
target_include_directories(pricer_lib PUBLIC include)

add_executable(convergence_test src/convergence_test.cpp)
target_link_libraries(convergence_test PRIVATE pricer_lib)

add_executable(run_hedge src/run_hedge.cpp)
target_link_libraries(run_hedge PRIVATE pricer_lib)

enable_testing()
add_subdirectory(tests)
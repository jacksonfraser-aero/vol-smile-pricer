add_executable(pricer_tests
    test_black_scholes.cpp
    test_binomial_tree.cpp
    test_delta_hedge.cpp
)
target_link_libraries(pricer_tests PRIVATE pricer_lib)
target_include_directories(pricer_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)

add_test(NAME pricer_tests COMMAND pricer_tests)


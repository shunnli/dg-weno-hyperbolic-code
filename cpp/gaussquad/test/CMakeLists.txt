add_executable(gaussquad_test gaussquad_test.cpp)
target_link_libraries(gaussquad_test PRIVATE gaussquad gtest_main)

include(GoogleTest)
gtest_discover_tests(gaussquad_test)

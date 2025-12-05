add_executable(gaussquad_test)
zero_get_files(tmp .)
target_sources(gaussquad_test PRIVATE ${tmp})
target_link_libraries(gaussquad_test PRIVATE gaussquad gtest_main)

include(GoogleTest)
gtest_discover_tests(gaussquad_test)

add_executable(example_fv_rk3_weno5_c)
target_sources(example_fv_rk3_weno5_c PRIVATE fv_rk3_weno5_c.cpp)
target_link_libraries(example_fv_rk3_weno5_c PRIVATE flux)
target_compile_definitions(example_fv_rk3_weno5_c PRIVATE OUTPUT_DIR="${EXAMPLE_OUTPUT_DIR}/FV-RK3-WENO5")
zero_check_target(example_fv_rk3_weno5_c)


add_executable(example_fv_rk3_weno5_v)
target_sources(example_fv_rk3_weno5_v PRIVATE fv_rk3_weno5_v.cpp)
target_link_libraries(example_fv_rk3_weno5_v PRIVATE flux)
target_compile_definitions(example_fv_rk3_weno5_v PRIVATE OUTPUT_DIR="${EXAMPLE_OUTPUT_DIR}/FV-RK3-WENO5")
zero_check_target(example_fv_rk3_weno5_v)

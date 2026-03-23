if(HAS_DEDUCING_THIS)
    add_executable(example_fv_godunov_n)
    target_sources(example_fv_godunov_n PRIVATE fv_godunov_n.cpp)
    target_link_libraries(example_fv_godunov_n PRIVATE flux)
    check_target(example_fv_godunov_n)
else()
    message(WARNING "Deducing this is not supported, skipping example_fv_godunov_n")
endif()


add_executable(example_fv_godunov_c)
target_sources(example_fv_godunov_c PRIVATE fv_godunov_c.cpp)
target_link_libraries(example_fv_godunov_c PRIVATE flux)
check_target(example_fv_godunov_c)

add_executable(example_fv_godunov_f)
target_sources(example_fv_godunov_f PRIVATE fv_godunov_f.cpp)
target_link_libraries(example_fv_godunov_f PRIVATE flux)
check_target(example_fv_godunov_f)

add_executable(example_fv_godunov_v)
target_sources(example_fv_godunov_v PRIVATE fv_godunov_v.cpp)
target_link_libraries(example_fv_godunov_v PRIVATE flux)
check_target(example_fv_godunov_v)

add_executable(example_fv_godunov_p)
target_sources(example_fv_godunov_p PRIVATE fv_godunov_p.cpp)
target_link_libraries(example_fv_godunov_p PRIVATE flux)
check_target(example_fv_godunov_p)

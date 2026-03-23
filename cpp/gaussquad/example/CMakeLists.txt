add_executable(quad_example quad_example.cpp)
target_link_libraries(quad_example PRIVATE gaussquad)

add_executable(quad3_example quad3_example.cpp)
target_link_libraries(quad3_example PRIVATE gaussquad)

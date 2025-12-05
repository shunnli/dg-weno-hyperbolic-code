cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(gaussquad VERSION 1.0)
include(cmake/zero.cmake)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

zero_setup()
zero_check()

add_library(gaussquad INTERFACE)
add_library(gaussquad::gaussquad ALIAS gaussquad)
target_include_directories(gaussquad INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(example)

    include(CTest)
    if(BUILD_TESTING)
        include(FetchContent)

        FetchContent_Declare(
            googletest
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
        )

        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        set(BUILD_GMOCK OFF CACHE BOOL "" FORCE) # Disable GMock
        set(INSTALL_GTEST OFF CACHE BOOL "" FORCE) # Disable install gtest

        FetchContent_MakeAvailable(googletest)

        add_subdirectory(test)
    endif()
endif()

install(TARGETS gaussquad EXPORT gaussquadTargets)
install(EXPORT gaussquadTargets
    FILE gaussquadConfig.cmake
    NAMESPACE gaussquad::
    DESTINATION lib/cmake/gaussquad)
install(DIRECTORY include/gaussquad DESTINATION include)

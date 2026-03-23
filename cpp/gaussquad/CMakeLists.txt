cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(gaussquad VERSION 1.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(PROJECT_BINARY_DIR STREQUAL PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "The binary directory cannot be the same as source directory.")
endif()

# create compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "" FORCE)

# c/c++ compile flags
if(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
    set(flags "-Wall -Wextra -Wfatal-errors -Wshadow -Wconversion -Wsign-conversion -Wuninitialized -pedantic -Wno-unused-parameter")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(flags "/W3 /WX /MP /utf-8 /permissive- /Zc:__cplusplus")
else()
    set(flags "")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags}")


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

        fetchcontent_declare(
            googletest
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
        )

        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        set(BUILD_GMOCK OFF CACHE BOOL "" FORCE) # Disable GMock
        set(INSTALL_GTEST OFF CACHE BOOL "" FORCE) # Disable install gtest

        fetchcontent_makeavailable(googletest)

        add_subdirectory(test)
    endif()
endif()

install(TARGETS gaussquad EXPORT gaussquadTargets)
install(EXPORT gaussquadTargets
    FILE gaussquadConfig.cmake
    NAMESPACE gaussquad::
    DESTINATION lib/cmake/gaussquad)
install(DIRECTORY include/gaussquad DESTINATION include)

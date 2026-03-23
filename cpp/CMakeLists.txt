cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(nnpde-cpp VERSION 1.0)

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


message(STATUS "------- <Check Information> -------")
message(STATUS ">> system = ${CMAKE_SYSTEM_NAME}")
message(STATUS ">> generator = ${CMAKE_GENERATOR}")
message(STATUS ">> build_type = ${CMAKE_BUILD_TYPE}")
message(STATUS ">> cxx_compiler_id = ${CMAKE_CXX_COMPILER_ID}(${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS ">> cxx_compiler = ${CMAKE_CXX_COMPILER}")
message(STATUS ">> cxx_flags = ${CMAKE_CXX_FLAGS}")
message(STATUS ">> cxx_flags_debug = ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS ">> cxx_flags_release = ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS ">> source_dir = ${PROJECT_SOURCE_DIR}")
message(STATUS ">> binary_dir = ${CMAKE_BINARY_DIR}")
message(STATUS ">> install_prefix = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "-----------------------------------")

include(cmake/check_target.cmake)


include(CheckCXXSourceCompiles)

# C++23 deducing this
check_cxx_source_compiles("
struct X {
    void f(this X& self) {}
};
int main() {
    X x; x.f();
    return 0;
}"    HAS_DEDUCING_THIS)

# std::format
check_cxx_source_compiles("
#include <format>
#include <string>
int main() {
    std::string s = std::format(\"Hello {}!\", 123);
    return 0;
}"    HAS_STD_FORMAT)

# std::expected
check_cxx_source_compiles("
#include <expected>
#include <string>
std::expected<int, std::string> func(bool ok) {
    if (ok) return 42;
    return std::unexpected(\"error\");
}
int main() {
    auto e = func(true);
    return e.value_or(-1);
}"    HAS_STD_EXPECTED)

message(STATUS "HAS_DEDUCING_THIS = ${HAS_DEDUCING_THIS}")
message(STATUS "HAS_STD_FORMAT    = ${HAS_STD_FORMAT}")
message(STATUS "HAS_STD_EXPECTED  = ${HAS_STD_EXPECTED}")

add_subdirectory(gaussquad)

add_library(flux INTERFACE)
target_include_directories(flux INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(flux INTERFACE gaussquad::gaussquad)
check_target(flux)

add_subdirectory(src)

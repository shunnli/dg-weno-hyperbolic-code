cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(nnpde-cpp VERSION 1.0)
include(cmake/zero.cmake)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

zero_setup()
zero_check()

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
zero_check_target(flux)

add_subdirectory(src)

#include <cmath>
#include <iomanip>
#include <iostream>

#include "gaussquad/gaussquad.hpp"

namespace {
// x^2 + y^2
double func1(double x, double y) { return x * x + y * y; }

// x * y
double func2(double x, double y) { return x * y; }

// sin(x) * cos(y)
double func3(double x, double y) { return std::sin(x) * std::cos(y); }
}  // namespace

using gaussquad::Quad3;
using gaussquad::Triangle;
using gaussquad::TriangleRule;

int main() {
    std::cout << std::setprecision(15);

    Triangle triangle1{
        .ax = 0.0, .ay = 0.0, .bx = 1.0, .by = 0.0, .cx = 0.0, .cy = 1.0};
    Triangle triangle2{
        .ax = 0.0, .ay = 0.0, .bx = 2.0, .by = 0.0, .cx = 0.0, .cy = 3.0};

    std::cout << "3 points Gauss Quadrature with func1 (x^2 + y^2): ";
    auto result = Quad3{TriangleRule::P3()}.integrate(func1, triangle1);
    std::cout << "\nInt(x^2 + y^2, Triangle) = " << result
              << " (expected 0.166667)\n";

    std::cout << "\n7 points Gauss Quadrature with func2 (x * y): ";
    result = Quad3{TriangleRule::P7()}.integrate(func2, triangle1);
    std::cout << "\nInt(x * y, Triangle) = " << result
              << " (expected 0.0416667)\n";

    std::cout << "\n7 points with func3 (sin(x) * cos(y)): ";
    result = Quad3{TriangleRule::P7()}.integrate(func3, triangle1);
    std::cout << "\nInt(sin(x) * cos(y), Triangle) = " << result
              << " (expected 0.150584)\n";

    std::cout << "\n12 points with func3 (sin(x) * cos(y)): ";
    result = Quad3{TriangleRule::P12()}.integrate(func3, triangle1);
    std::cout << "\nInt(sin(x) * cos(y), Triangle) = " << result
              << " (expected 0.150584)\n";

    std::cout << "\nTriangle (0, 0), (2, 0), (0, 3) with func1 (x^2 + y^2): ";
    result = Quad3{TriangleRule::P12()}.integrate(func1, triangle2);
    std::cout << "\nInt(x^2 + y^2, Custom Triangle) = " << result
              << " (expected 6.5)\n";

    return 0;
}

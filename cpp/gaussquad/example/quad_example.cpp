#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numbers>

#include "gaussquad/gaussquad.hpp"

using gaussquad::gausslegendre;
using gaussquad::gausslobatto;
using gaussquad::Quad;

double func(double x) { return x * x; }

int main() {
    std::cout << std::setprecision(15);

    std::cout << "Gauss-Legendre 3 points: ";
    double result =
        Quad(gausslegendre(3)).integrate(std::function<double(double)>(func));
    std::cout << "NInt(x^2,{x,-1,1}) = " << result << "\n";

    std::cout << "\nGauss-Lobatto 5 points: ";
    result =
        Quad(gausslobatto(5)).integrate(std::function<double(double)>(func));
    std::cout << "NInt(x^2,{x,-1,1}) = " << result << "\n";

    std::cout << "\nGauss-Legendre 5 points: ";
    auto func1 = [](double x) { return std::sin(x); };
    result = Quad(gausslegendre(5)).integrate(func1);
    std::cout << "NInt(sin(x),{x,-1,1}) = " << result << "\n";

    // Mapping [-1,1] -> [0,pi]
    auto pi_f  = [](double x) { return std::numbers::pi / 2 * (x + 1); };
    auto pi_df = [](double)   { return std::numbers::pi / 2; };

    std::cout << "\nGauss-Legendre 5 points: ";
    result = Quad(gausslegendre(5))
                 .transform(pi_f, pi_df)
                 .integrate(std::function<double(double)>(func1));
    std::cout << "NInt(sin(x),{x,0,pi}) = " << result << "\n";

    std::cout << "\nGauss-Legendre 13 points: ";
    result = Quad(gausslegendre(13)).transform(pi_f, pi_df).integrate(func1);
    std::cout << "NInt(sin(x),{x,0,pi}) = " << result << "\n";

    std::cout << "\nGauss-Lobatto 11 points: ";
    result = Quad(gausslobatto(11)).transform(pi_f, pi_df).integrate(func1);
    std::cout << "NInt(sin(x),{x,0,pi}) = " << result << "\n";

    return 0;
}

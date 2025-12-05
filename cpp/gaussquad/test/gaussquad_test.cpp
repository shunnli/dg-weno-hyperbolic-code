#include "gtest/gtest.h"

#include "gaussquad/gaussquad.hpp"

#include <cmath>

using gaussquad::gausslegendre;
using gaussquad::gausslobatto;

double reference_integral(unsigned degree, double xl, double xr) {
    return (std::pow(xr, degree + 1) - std::pow(xl, degree + 1))
           / static_cast<double>(degree + 1);
}

TEST(GaussQuadratureTest, GaussLegendreAccuracy) {
    for (unsigned order = 2; order <= 18; ++order) {
        auto [points, weights] = gausslegendre(order);

        for (unsigned degree = 0; degree <= 2 * order - 1; ++degree) {
            double integral = 0.0;
            for (unsigned i = 0; i < order; ++i) {
                integral += weights[i] * std::pow(points[i], degree);
            }

            double exact = reference_integral(degree, -1.0, 1.0);
            EXPECT_NEAR(integral, exact, 1e-12)
                << "Failed for polynomial degree " << degree;
        }
    }
}

TEST(GaussQuadratureTest, GaussLobattoAccuracy) {
    for (unsigned order = 2; order <= 18; ++order) {
        auto [points, weights] = gausslobatto(order);

        for (unsigned degree = 0; degree <= 2 * order - 3; ++degree) {
            double integral = 0.0;
            for (unsigned i = 0; i < order; ++i) {
                integral += weights[i] * std::pow(points[i], degree);
            }

            double exact = reference_integral(degree, -1.0, 1.0);
            EXPECT_NEAR(integral, exact, 1e-12)
                << "Failed for polynomial degree " << degree;
        }
    }
}

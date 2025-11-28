#pragma once

#include <functional>
#include <vector>

#include "burgers_exact.hpp"
#include "constants.hpp"

using namespace flux;      // NOLINT
using flux::constant::pi;  // NOLINT

struct Config {
    double xl;
    double xr;
    double tend;
    size_t gauss_k;
    std::vector<size_t> nlist;

    std::function<double(double)> init;
    std::function<double(double, double)> exact;
};

inline auto plot_config() {
    return Config{
        .xl = -pi,
        .xr = pi,
        .tend = 1.5,
        .gauss_k = 5,
        .nlist = {20, 80},
        .init =
            [](double x) {
                return BurgersExact(0.5, 1.0, 1.0, 0, 1e-10).eval(x, 0);
            },
        .exact =
            [](double x, double t) {
                return BurgersExact(0.5, 1.0, 1.0, 0, 1e-10).eval(x, t);
            },
    };
}

inline auto order_test_config() {
    return Config{
        .xl = -pi,
        .xr = pi,
        .tend = 0.5,
        .gauss_k = 5,
        .nlist = {10, 20, 40, 80, 160, 320, 640},
        .init =
            [](double x) {
                return BurgersExact(0.5, 1.0, 1.0, 0, 1e-10).eval(x, 0);
            },
        .exact =
            [](double x, double t) {
                return BurgersExact(0.5, 1.0, 1.0, 0, 1e-10)
                    .eval_with_check(x, t);
            },
    };
}

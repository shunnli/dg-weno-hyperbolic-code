#pragma once

#include <vector>

#include "period_index.hpp"

namespace flux {
inline void weno5(std::vector<double> u, std::vector<double> &res_ul,
                  std::vector<double> &res_ur) {
    // linear weight
    constexpr double d_l0 = 3.0 / 10;
    constexpr double d_l1 = 3.0 / 5;
    constexpr double d_l2 = 1.0 / 10;
    constexpr double d_r0 = 1.0 / 10;
    constexpr double d_r1 = 3.0 / 5;
    constexpr double d_r2 = 3.0 / 10;

    // Important, avoid denominator being 0 and not too small
    constexpr double weno_ep = 1e-6;

    const auto cb = [](double u0, double u1, double u2, double c0, double c1,
                       double c2) { return c0 * u0 + c1 * u1 + c2 * u2; };
    const auto cb2 = [](double u0, double u1, double c0, double c1) {
        return c0 * u0 * u0 + c1 * u1 * u1;
    };

    size_t n = u.size();
    res_ul = std::vector<double>(n);
    res_ur = std::vector<double>(n);
    for (size_t i = 0; i < n; i++) {
        auto idx = PeriodIndex(n, i);

        // smooth indicator
        double b0 = cb2(cb(u[idx.l(2)], u[idx.l()], u[idx.c()], 1, -2, 1),
                        cb(u[idx.l(2)], u[idx.l()], u[idx.c()], 1, -4, 3),  //
                        13.0 / 12, 1.0 / 4);
        double b1 = cb2(cb(u[idx.l()], u[idx.c()], u[idx.r()], 1, -2, 1),
                        cb(u[idx.l()], u[idx.c()], u[idx.r()], 1, 0, -1),  //
                        13.0 / 12, 1.0 / 4);
        double b2 = cb2(cb(u[idx.c()], u[idx.r()], u[idx.r(2)], 1, -2, 1),
                        cb(u[idx.c()], u[idx.r()], u[idx.r(2)], 3, -4, 1),  //
                        13.0 / 12, 1.0 / 4);

        // Nonlinear weight
        double a_l0 = d_l0 / ((b0 + weno_ep) * (b0 + weno_ep));
        double a_l1 = d_l1 / ((b1 + weno_ep) * (b1 + weno_ep));
        double a_l2 = d_l2 / ((b2 + weno_ep) * (b2 + weno_ep));
        double a_r0 = d_r0 / ((b0 + weno_ep) * (b0 + weno_ep));
        double a_r1 = d_r1 / ((b1 + weno_ep) * (b1 + weno_ep));
        double a_r2 = d_r2 / ((b2 + weno_ep) * (b2 + weno_ep));

        // Normalized nonlinear weight
        double a_l_sum = a_l0 + a_l1 + a_l2;
        double a_r_sum = a_r0 + a_r1 + a_r2;
        double w_l0 = a_l0 / a_l_sum;
        double w_l1 = a_l1 / a_l_sum;
        double w_l2 = a_l2 / a_l_sum;
        double w_r0 = a_r0 / a_r_sum;
        double w_r1 = a_r1 / a_r_sum;
        double w_r2 = a_r2 / a_r_sum;

        double u_l0 = cb(u[idx.l(2)], u[idx.l()], u[idx.c()],  //
                         -1.0 / 6, 5.0 / 6, 1.0 / 3);
        double u_l1 = cb(u[idx.l()], u[idx.c()], u[idx.r()],  //
                         1.0 / 3, 5.0 / 6, -1.0 / 6);
        double u_l2 = cb(u[idx.c()], u[idx.r()], u[idx.r(2)],  //
                         11.0 / 6, -7.0 / 6, 1.0 / 3);

        double u_r0 = cb(u[idx.l(2)], u[idx.l()], u[idx.c()],  //
                         1.0 / 3, -7.0 / 6, 11.0 / 6);
        double u_r1 = cb(u[idx.l()], u[idx.c()], u[idx.r()],  //
                         -1.0 / 6, 5.0 / 6, 1.0 / 3);
        double u_r2 = cb(u[idx.c()], u[idx.r()], u[idx.r(2)],  //
                         1.0 / 3, 5.0 / 6, -1.0 / 6);

        res_ul[idx.c()] = w_l0 * u_l0 + w_l1 * u_l1 + w_l2 * u_l2;
        res_ur[idx.c()] = w_r0 * u_r0 + w_r1 * u_r1 + w_r2 * u_r2;
    }
    return;
}
}  // namespace flux

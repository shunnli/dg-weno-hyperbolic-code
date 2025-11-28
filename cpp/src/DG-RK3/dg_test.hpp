#include "config.hpp"
#include "legendre_polys.hpp"
#include "linespace.hpp"

#include "error_and_order.hpp"
#include "export_to_file.hpp"

#include "solver/preset.hpp"

#include "gaussquad/gausslegendre.hpp"

using namespace flux;  // NOLINT

using P = LegendrePolys;
using Px = LegendrePolysDx;

inline std::vector<double>
DG_projection(const std::function<double(double)> &u0,
              const std::vector<double> &x, double dx, size_t DG_k,
              size_t gauss_k) {
    size_t cell_num = x.size();
    std::vector<double> uh(cell_num * (DG_k + 1));

    auto [gauss_points, gauss_weights] =
        gaussquad::gausslegendre(static_cast<unsigned>(gauss_k));

    for (size_t i = 0; i < cell_num; i++) {
        for (size_t j = 0; j <= DG_k; j++) {
            double tmp_sum = 0;
            for (size_t gauss_i = 0; gauss_i < gauss_k; gauss_i++) {
                double tmp1 = u0(x[i] + gauss_points[gauss_i] * dx / 2);
                double tmp2 = P::eval(j, gauss_points[gauss_i]);
                tmp_sum += gauss_weights[gauss_i] * tmp1 * tmp2;
            }
            uh[i * (DG_k + 1) + j] =
                tmp_sum * static_cast<double>(2 * j + 1) / 2;
        }
    }

    return uh;
}

inline auto DG_error(const std::vector<double> &uh,
                     const std::function<double(double)> &uexact,
                     const std::vector<double> &x, double dx, size_t DG_k,
                     size_t gauss_k) {
    size_t cell_num = x.size();

    auto [gauss_points, gauss_weights] =
        gaussquad::gausslegendre(static_cast<unsigned>(gauss_k));

    double error_linf = 0;
    double error_l1 = 0;
    double error_l2 = 0;
    for (size_t i = 0; i < cell_num; ++i) {
        for (size_t gauss_i = 0; gauss_i < gauss_k; gauss_i++) {
            double uh_value =
                evals<P>(uh, gauss_points[gauss_i], i * (DG_k + 1), DG_k + 1);
            double uexact_value = uexact(x[i] + dx / 2 * gauss_points[gauss_i]);

            double tmp = std::abs(uh_value - uexact_value);

            if (tmp > error_linf) error_linf = tmp;

            error_l1 += gauss_weights[gauss_i] * tmp * dx / 2;

            error_l2 += gauss_weights[gauss_i] * tmp * tmp * dx / 2;
        }
    }

    error_l2 = std::sqrt(error_l2);

    return std::make_tuple(error_l1, error_l2, error_linf);
}

template <typename SolverType>
void DG_plot_test(Config cfg, SolverType solver, size_t DG_k,
                  const std::vector<const char *> &filelist) {
    double dx = 0;

    for (size_t i = 0; i < cfg.nlist.size(); i++) {
        size_t n = cfg.nlist[i];
        auto x = linespace_mid(cfg.xl, cfg.xr, n, dx);

        // L2 Projection
        auto uh = DG_projection(cfg.init, x, dx, DG_k, cfg.gauss_k);

        auto ex = Mesh1d{dx};
        uh = solver.run(Vec{uh}, ex, 0, cfg.tend).value().data;

        // midpoint value
        auto uh_data = std::vector<double>(n);
        auto u_data = std::vector<double>(n);
        for (size_t j = 0; j < n; j++) {
            uh_data[j] = evals<P>(uh, 0, j * (DG_k + 1), DG_k + 1);
            u_data[j] = cfg.exact(x[j], cfg.tend);
        }

        export_to_file(filelist[i], x, u_data, uh_data, ',');
    }
}

template <typename SolverType>
void DG_order_test(Config cfg, SolverType solver, size_t DG_k,
                   const char *filename) {
    double dx = 0;
    size_t gauss_k = cfg.gauss_k;

    auto [gauss_points, _] = gaussquad::gausslegendre(static_cast<unsigned>(gauss_k));

    auto &nlist = cfg.nlist;

    auto error_l1 = std::vector<double>(nlist.size());
    auto error_l2 = std::vector<double>(nlist.size());
    auto error_linf = std::vector<double>(nlist.size());

    for (size_t i = 0; i < nlist.size(); i++) {
        size_t n = nlist[i];
        auto x = linespace_mid(cfg.xl, cfg.xr, n, dx);

        // L2 Projection
        auto uh = DG_projection(cfg.init, x, dx, DG_k, gauss_k);

        auto ex = Mesh1d{dx};
        uh = solver.run(Vec{uh}, ex, 0, cfg.tend).value().data;

        auto errs = DG_error(
            uh, [cfg](double s) { return cfg.exact(s, cfg.tend); }, x, dx, DG_k,
            gauss_k);

        error_l1[i] = std::get<0>(errs);
        error_l2[i] = std::get<1>(errs);
        error_linf[i] = std::get<2>(errs);
    }

    auto order_l1 = order(error_l1, nlist);
    auto order_l2 = order(error_l2, nlist);
    auto order_linf = order(error_linf, nlist);

    print_error_table(std::cout, nlist, error_l1, error_l2, error_linf,
                      order_l1, order_l2, order_linf, ' ');
    print_error_table_to_file(filename, nlist, error_l1, error_l2, error_linf,
                              order_l1, order_l2, order_linf, '&');
}

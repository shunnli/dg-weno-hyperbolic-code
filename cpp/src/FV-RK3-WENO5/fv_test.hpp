#include "config.hpp"
#include "linespace.hpp"

#include "error_and_order.hpp"
#include "export_to_file.hpp"

#include "solver/preset.hpp"

#include "gaussquad/gausslegendre.hpp"
#include "gaussquad/quad.hpp"

using namespace flux;  // NOLINT

template <typename SolverType>
void FV_plot_test(Config cfg, SolverType solver,
                  const std::vector<std::filesystem::path> &filelist) {
    double dx = 0;

    auto exact = [=](double x) { return cfg.exact(x, cfg.tend); };

    auto g = gaussquad::Quad(
        gaussquad::gausslegendre(static_cast<unsigned>(cfg.gauss_k)));

    for (size_t i = 0; i < cfg.nlist.size(); i++) {
        size_t n = cfg.nlist[i];
        auto x = linespace_mid(cfg.xl, cfg.xr, n, dx);
        auto uh = std::vector<double>(n);

        for (size_t j = 0; j < n; j++) {
            double tmp = g.integrate([=](double s) {
                return cfg.init(x[j] + s * dx / 2);
            }) * dx / 2;
            uh[j] = tmp / dx;
        }

        auto ex = Mesh1d{dx};
        uh = solver.run(Vec{uh}, ex, 0, cfg.tend).value().data;

        auto u = std::vector<double>(n);
        for (size_t j = 0; j < n; j++) {
            double tmp = g.integrate([=](double s) {
                return exact(x[j] + s * dx / 2);
            }) * dx / 2;
            u[j] = tmp / dx;
        }

        export_to_file(filelist[i], x, u, uh, ',');
    }
}

template <typename SolverType>
void FV_order_test(Config cfg, SolverType solver,
                   const std::filesystem::path &filename) {
    double dx = 0;

    auto exact = [=](double x) { return cfg.exact(x, cfg.tend); };

    auto g = gaussquad::Quad(
        gaussquad::gausslegendre(static_cast<unsigned>(cfg.gauss_k)));

    auto error_l1 = std::vector<double>(cfg.nlist.size());
    auto error_l2 = std::vector<double>(cfg.nlist.size());
    auto error_linf = std::vector<double>(cfg.nlist.size());

    for (size_t i = 0; i < cfg.nlist.size(); i++) {
        size_t n = cfg.nlist[i];
        auto x = linespace_mid(cfg.xl, cfg.xr, n, dx);
        auto uh = std::vector<double>(n);

        for (size_t j = 0; j < n; j++) {
            double tmp = g.integrate([=](double s) {
                return cfg.init(x[j] + s * dx / 2);
            }) * dx / 2;
            uh[j] = tmp / dx;
        }

        auto ex = Mesh1d{dx};
        uh = solver.run(Vec{uh}, ex, 0, cfg.tend).value().data;

        auto u = std::vector<double>(n);
        for (size_t j = 0; j < n; j++) {
            double tmp = g.integrate([=](double s) {
                return exact(x[j] + s * dx / 2);
            }) * dx / 2;
            u[j] = tmp / dx;
        }

        error_l1[i] = error(uh, u, dx, ErrorType::L1);
        error_l2[i] = error(uh, u, dx, ErrorType::L2);
        error_linf[i] = error(uh, u, dx, ErrorType::Linf);
    }

    auto order_l1 = order(error_l1, cfg.nlist);
    auto order_l2 = order(error_l2, cfg.nlist);
    auto order_linf = order(error_linf, cfg.nlist);

    print_error_table(std::cout, cfg.nlist, error_l1, error_l2, error_linf,
                      order_l1, order_l2, order_linf, ' ');
    print_error_table_to_file(filename, cfg.nlist, error_l1, error_l2,
                              error_linf, order_l1, order_l2, order_linf, '&');
}

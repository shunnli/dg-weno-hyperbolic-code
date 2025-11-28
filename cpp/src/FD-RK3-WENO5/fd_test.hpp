#include "config.hpp"
#include "linespace.hpp"

#include "error_and_order.hpp"
#include "export_to_file.hpp"

#include "solver/preset.hpp"

using namespace flux;  // NOLINT

template <typename SolverType>
void FD_plot_test(Config cfg, SolverType solver,
                  const std::vector<const char *> &filelist) {
    double dx = 0;

    for (size_t i = 0; i < cfg.nlist.size(); i++) {
        size_t n = cfg.nlist[i];
        auto x = linespace_mid(cfg.xl, cfg.xr, n, dx);
        auto uh = std::vector<double>(n);

        for (size_t j = 0; j < n; j++) { uh[j] = cfg.init(x[j]); }

        auto ex = Mesh1d{dx};
        uh = solver.run(Vec{uh}, ex, 0, cfg.tend).value().data;

        auto u = std::vector<double>(n);
        for (size_t j = 0; j < n; j++) { u[j] = cfg.exact(x[j], cfg.tend); }

        export_to_file(filelist[i], x, u, uh, ',');
    }
}

template <typename SolverType>
void FD_order_test(Config cfg, SolverType solver, const char *filename) {
    double dx = 0;

    auto error_l1 = std::vector<double>(cfg.nlist.size());
    auto error_l2 = std::vector<double>(cfg.nlist.size());
    auto error_linf = std::vector<double>(cfg.nlist.size());

    for (size_t i = 0; i < cfg.nlist.size(); i++) {
        size_t n = cfg.nlist[i];
        auto x = linespace_mid(cfg.xl, cfg.xr, n, dx);
        auto uh = std::vector<double>(n);

        for (size_t j = 0; j < n; j++) { uh[j] = cfg.init(x[j]); }

        auto ex = Mesh1d{dx};
        uh = solver.run(Vec{uh}, ex, 0, cfg.tend).value().data;

        auto u = std::vector<double>(n);
        for (size_t j = 0; j < n; j++) { u[j] = cfg.exact(x[j], cfg.tend); }

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

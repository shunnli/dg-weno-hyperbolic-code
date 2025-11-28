#include "fv_test.hpp"
#include "period_index.hpp"
#include "solver/solver_stdfunc.hpp"

using namespace flux;  // NOLINT
using flux::solver_stdfunc::Solver;
using flux::solver_stdfunc::UpdaterFactory;

auto FV_godunov_solver() {
    auto solver = Solver<Vec, Mesh1d>{};

    auto get_dt = [](const Vec &var, Mesh1d &ex, double t) {
        const auto &u = var.data;
        double df_max = 0;
        for (const auto ui : u) {
            double tmp = std::abs(ui);  // df(u) = u
            if (tmp > df_max) df_max = tmp;
        }
        return 0.5 * ex.dx / df_max;
    };

    auto fhat_godunov = [](double ul, double ur) {
        if (ul <= ur) {  // min
            if (ul * ur > 0) { return std::min(ul * ul / 2, ur * ur / 2); }
            return 0.0;
        }
        return std::max(ul * ul / 2, ur * ur / 2);  // max
    };

    auto op_L = [=](const Vec &var, Mesh1d &ex, double t) {
        const auto &u = var.data;
        size_t n = u.size();
        auto L = std::vector<double>(n);
        for (size_t i = 0; i < n; i++) {
            auto idx = PeriodIndex(n, i);
            double fhat_l = fhat_godunov(u[idx.l()], u[idx.c()]);
            double fhat_r = fhat_godunov(u[idx.c()], u[idx.r()]);
            L[i] = (fhat_l - fhat_r) / ex.dx;
        }
        return Vec{L};
    };

    auto euler_updater =
        UpdaterFactory<Vec, Mesh1d>::get_euler_updater(op_L, get_dt, {}, {});

    return solver.set_update(euler_updater);
}

int main() {
    auto solver = FV_godunov_solver();
    FV_order_test(order_test_config(), solver, OUTPUT_DIR "/order_f.csv");
    FV_plot_test(plot_config(), solver,
                 {OUTPUT_DIR "/plot_1_f.csv", OUTPUT_DIR "/plot_2_f.csv"});

    return 0;
}

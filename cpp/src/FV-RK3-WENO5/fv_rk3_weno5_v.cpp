#include <algorithm>

#include "fv_test.hpp"
#include "period_index.hpp"
#include "solver/solver_virtual.hpp"
#include "weno5.hpp"

#include "output_manager.hpp"

using namespace flux;                   // NOLINT
using flux::solver_virtual::RK3Solver;  // NOLINT

class FVWENO5Solver : public RK3Solver<Vec, Mesh1d> {
public:
    double get_dt(const Vec &var, Mesh1d &ex, double t) const override {
        double df_max = 0;
        for (const auto ui : var.data) {
            double tmp = std::abs(ui);  // df(u) = u
            df_max = std::max(tmp, df_max);
        }
        return std::pow(ex.dx, 5.0 / 3) / (2 * df_max);
    }

    Vec op_L(const Vec &var, Mesh1d &ex, double t) const override {
        const auto &u = var.data;
        size_t n = u.size();
        auto L = std::vector<double>(n);
        auto ul_p = std::vector<double>(n);
        auto ur_m = std::vector<double>(n);

        weno5(u, ul_p, ur_m);  // WENO

        for (size_t i = 0; i < n; i++) {
            auto idx = PeriodIndex(n, i);
            double fhat_l = fhat_LF(ur_m[idx.l()], ul_p[idx.c()]);
            double fhat_r = fhat_LF(ur_m[idx.c()], ul_p[idx.r()]);
            L[i] = (fhat_l - fhat_r) / ex.dx;
        }
        return Vec{L};
    }

    static double fhat_LF(double ul, double ur) {
        double c = std::max(std::abs(ul), std::abs(ur));

        double tmp1 = 0.5 * (ul * ul / 2 + ur * ur / 2);
        double tmp2 = 0.5 * c * (ur - ul);
        return tmp1 - tmp2;
    };
};

int main(int argc, char **argv) {
    auto output_path = OutputManager::parse_output(argc, argv);
    const auto out =
        OutputManager(output_path.value_or("outputs"), "FV-RK3-WENO5");

    auto solver = FVWENO5Solver{};
    FV_order_test(order_test_config(), solver, out / "order_v.csv");
    FV_plot_test(plot_config(), solver,
                 {out / "plot_1_v.csv", out / "plot_2_v.csv"});

    return 0;
}

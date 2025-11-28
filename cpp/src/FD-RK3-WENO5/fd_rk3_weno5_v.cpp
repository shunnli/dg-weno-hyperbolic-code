#include "fd_test.hpp"
#include "period_index.hpp"
#include "solver/solver_virtual.hpp"
#include "weno5.hpp"

using namespace flux;  // NOLINT
using flux::solver_virtual::RK3Solver;

class FDWENO5Solver : public RK3Solver<Vec, Mesh1d> {
public:
    double get_dt(const Vec &var, Mesh1d &ex, double t) const override {
        double df_max = 0;
        for (const auto ui : var.data) {
            double tmp = std::abs(ui);  // df(u) = u
            if (tmp > df_max) df_max = tmp;
        }
        return std::pow(ex.dx, 5.0 / 3) / (2 * df_max);
    }

    Vec op_L(const Vec &var, Mesh1d &ex, double t) const override {
        const auto &u = var.data;
        size_t n = u.size();
        auto L = std::vector<double>(n);

        double lf_c{0};  // gobal c
        for (size_t i = 0; i < n; i++) {
            double lf_c_tmp = std::abs(u[i]);
            lf_c = (lf_c_tmp > lf_c) ? lf_c_tmp : lf_c;
        }

        // split
        auto fplus = [lf_c](double v) { return 0.5 * (v * v / 2 + lf_c * v); };
        auto fminus = [lf_c](double v) { return 0.5 * (v * v / 2 - lf_c * v); };
        auto fu_plus = std::vector<double>(n);
        auto fu_minus = std::vector<double>(n);

        for (size_t i = 0; i < n; i++) {
            fu_plus[i] = fplus(u[i]);
            fu_minus[i] = fminus(u[i]);
        }

        auto fplus_r = std::vector<double>(n);
        auto fplus_l_useless = std::vector<double>(n);   // useless
        auto fminus_r_useless = std::vector<double>(n);  // useless
        auto fminus_l = std::vector<double>(n);

        weno5(fu_plus, fplus_l_useless, fplus_r);
        weno5(fu_minus, fminus_l, fminus_r_useless);

        for (size_t i = 0; i < n; i++) {
            auto idx = PeriodIndex(n, i);
            double fhat_l = fplus_r[idx.l()] + fminus_l[idx.c()];
            double fhat_r = fplus_r[idx.c()] + fminus_l[idx.r()];
            L[i] = (fhat_l - fhat_r) / ex.dx;
        }
        return Vec{L};
    }
};

int main() {
    auto solver = FDWENO5Solver{};
    FD_order_test(order_test_config(), solver, OUTPUT_DIR "/order_v.csv");
    FD_plot_test(plot_config(), solver,
                 {OUTPUT_DIR "/plot_1_v.csv", OUTPUT_DIR "/plot_2_v.csv"});

    return 0;
}

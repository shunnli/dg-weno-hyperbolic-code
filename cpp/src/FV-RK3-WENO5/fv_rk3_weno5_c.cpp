#include "fv_test.hpp"
#include "period_index.hpp"
#include "solver/solver_crtp.hpp"
#include "weno5.hpp"

using namespace flux;                // NOLINT
using flux::solver_crtp::RK3Solver;  // NOLINT

class FVWENO5Solver : public RK3Solver<Vec, Mesh1d, FVWENO5Solver> {
public:
    static double get_dt(const Vec &var, Mesh1d &ex, double t) {
        double df_max = 0;
        for (const auto ui : var.data) {
            double tmp = std::abs(ui);  // df(u) = u
            if (tmp > df_max) df_max = tmp;
        }
        return std::pow(ex.dx, 5.0 / 3) / (2 * df_max);
    }

    static Vec op_L(const Vec &var, Mesh1d &ex, double t) {
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

int main() {
    auto solver = FVWENO5Solver{};
    FV_order_test(order_test_config(), solver, OUTPUT_DIR "/order_c.csv");
    FV_plot_test(plot_config(), solver,
                 {OUTPUT_DIR "/plot_1_c.csv", OUTPUT_DIR "/plot_2_c.csv"});

    return 0;
}

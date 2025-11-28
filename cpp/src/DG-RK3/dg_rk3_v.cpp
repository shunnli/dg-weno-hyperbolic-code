#include "dg_test.hpp"

#include "legendre_polys.hpp"
#include "limiter.hpp"
#include "period_index.hpp"
#include "solver/solver_virtual.hpp"

#include "gaussquad/gausslegendre.hpp"

using namespace flux;  // NOLINT
using flux::solver_virtual::RK3Solver;

using P = LegendrePolys;
using Px = LegendrePolysDx;

class DGSolver : public RK3Solver<Vec, Mesh1d> {
public:
    size_t m_DG_k;
    size_t m_gauss_k;

    DGSolver(size_t DG_k, size_t gauss_k) : m_DG_k(DG_k), m_gauss_k(gauss_k) {}

    double get_dt(const Vec &var, Mesh1d &ex, double t) const override {
        const auto &u = var.data;

        double df_max = 0;
        size_t cell_num = u.size() / (m_DG_k + 1);

        for (size_t i = 0; i < cell_num; i++) {
            double tmp = std::abs(evals<P>(u, 0, i * (m_DG_k + 1), m_DG_k + 1));
            if (tmp > df_max) df_max = tmp;
        }

        auto coeff = static_cast<double>(2 * m_DG_k + 1);  // DG CFL

        if (m_DG_k > 2) {
            return pow(ex.dx, static_cast<double>(m_DG_k + 1) / 3)
                   / (coeff * df_max);
        }
        return ex.dx / (coeff * df_max);
    }

    Vec op_L(const Vec &var, Mesh1d &ex, double t) const override {
        const auto &u = var.data;
        size_t cell_num = u.size() / (m_DG_k + 1);

        auto ul = std::vector<double>(cell_num);
        auto uc = std::vector<double>(cell_num);
        auto ur = std::vector<double>(cell_num);
        for (size_t i = 0; i < cell_num; i++) {
            ul[i] = evals<P>(u, -1, i * (m_DG_k + 1), m_DG_k + 1);
            uc[i] = evals<P>(u, 0, i * (m_DG_k + 1), m_DG_k + 1);
            ur[i] = evals<P>(u, 1, i * (m_DG_k + 1), m_DG_k + 1);
        }

        auto fhat_l = std::vector<double>(cell_num);
        auto fhat_r = std::vector<double>(cell_num);

        for (size_t i = 0; i < cell_num; i++) {
            auto idx = PeriodIndex(cell_num, i);

            fhat_l[idx.c()] = fhat_LF(ur[idx.l()], ul[idx.c()]);
            fhat_r[idx.c()] = fhat_LF(ur[idx.c()], ul[idx.r()]);
        }

        auto [gauss_points, gauss_weights] =
            gaussquad::gausslegendre(static_cast<unsigned>(m_gauss_k));

        auto L = std::vector<double>(u.size());
        for (size_t i = 0; i < cell_num; i++) {
            for (size_t j = 0; j <= m_DG_k; j++) {
                double tmp_sum = 0;
                for (size_t gauss_i = 0; gauss_i < m_gauss_k; gauss_i++) {
                    double tmp1 = evals<P>(u, gauss_points[gauss_i],
                                           i * (m_DG_k + 1), m_DG_k + 1);
                    double tmp2 = tmp1 * tmp1 / 2;
                    double tmp3 =
                        Px::eval(j, gauss_points[gauss_i]) * (2 / ex.dx);
                    tmp_sum += gauss_weights[gauss_i] * tmp2 * tmp3;
                }
                double Fu = tmp_sum * (ex.dx / 2);
                double bl = fhat_l[i] * P::eval(j, -1);
                double br = fhat_r[i] * P::eval(j, 1);

                double inner_inv = static_cast<double>(2 * j + 1) / ex.dx;

                L[i * (m_DG_k + 1) + j] = inner_inv * (Fu - br + bl);
            }
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

class DGSolverWithLimiter : public DGSolver {
public:
    double m_tvb_M;

    DGSolverWithLimiter(size_t DG_k, size_t gauss_k, double tvb_M)
        : DGSolver(DG_k, gauss_k), m_tvb_M(tvb_M) {}

    Vec post_process_rk_stage(const Vec &var, Mesh1d &ex,
                              double t) const override {
        const auto &u = var.data;
        size_t cell_num = u.size() / (m_DG_k + 1);

        auto ul = std::vector<double>(cell_num);
        auto u_mean = std::vector<double>(cell_num);
        auto ur = std::vector<double>(cell_num);
        for (size_t i = 0; i < cell_num; i++) {
            ul[i] = evals<P>(u, -1, i * (m_DG_k + 1), m_DG_k + 1);
            u_mean[i] = u[i * (m_DG_k + 1)];
            ur[i] = evals<P>(u, 1, i * (m_DG_k + 1), m_DG_k + 1);
        }

        auto limiter = Limiter{m_tvb_M * ex.dx * ex.dx};  // add limiter

        auto u2 = std::vector<double>(u);
        for (size_t i = 0; i < cell_num; i++) {
            auto idx = PeriodIndex(cell_num, i);

            double ret_ul = ul[idx.c()];
            double ret_ur = ur[idx.c()];
            limiter.minmod(ret_ul, ret_ur, u_mean[idx.l()], u_mean[idx.c()],
                           u_mean[idx.r()]);

            Limiter::DG_recover(u2, i * (m_DG_k + 1), m_DG_k, u_mean[idx.c()],
                                ret_ul, ret_ur);
        }

        return Vec{u2};
    }
};

int main() {
    size_t DG_k = 2;
    size_t gauss_k = 7;
    auto cig_o = order_test_config();
    cig_o.gauss_k = gauss_k;
    auto cfg_p = plot_config();
    cfg_p.gauss_k = gauss_k;

    // no limiter

    auto solver1 = DGSolver{DG_k, gauss_k};
    DG_order_test(cig_o, solver1, DG_k, OUTPUT_DIR "/order_1_v.csv");
    DG_plot_test(cfg_p, solver1, DG_k,
                 {OUTPUT_DIR "/plot_11_v.csv", OUTPUT_DIR "/plot_12_v.csv"});

    // with limiter

    auto solver2 = DGSolverWithLimiter{DG_k, gauss_k, 0};
    DG_order_test(cig_o, solver2, DG_k, OUTPUT_DIR "/order_2_v.csv");
    DG_plot_test(cfg_p, solver2, DG_k,
                 {OUTPUT_DIR "/plot_21_v.csv", OUTPUT_DIR "/plot_22_v.csv"});

    auto solver3 = DGSolverWithLimiter{DG_k, gauss_k, 1.0};
    DG_order_test(cig_o, solver3, DG_k, OUTPUT_DIR "/order_3_v.csv");
    DG_plot_test(cfg_p, solver3, DG_k,
                 {OUTPUT_DIR "/plot_31_v.csv", OUTPUT_DIR "/plot_32_v.csv"});

    return 0;
}

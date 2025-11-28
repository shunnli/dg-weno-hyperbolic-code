#pragma once

#include <cmath>
#include <vector>

namespace flux {
class Limiter {
public:
    explicit Limiter(double tvb_M) : m_tvb_M(tvb_M) {}

    void minmod(double &ret_uleft_p, double &ret_uright_m, double uleft_mean,
                double u_mean, double uright_mean) const {
        bool flag{false};

        auto tmp1 = minmod_kernel(ret_uright_m - u_mean, uright_mean - u_mean,
                                  u_mean - uleft_mean, flag);
        auto tmp2 = minmod_kernel(u_mean - ret_uleft_p, uright_mean - u_mean,
                                  u_mean - uleft_mean, flag);

        ret_uright_m = u_mean + tmp1;
        ret_uleft_p = u_mean - tmp2;
    }

    double minmod_kernel(double a1, double a2, double a3, bool &flag) const {
        if (std::fabs(a1) < m_tvb_M) return a1;

        if (a1 > 0 && a2 > 0 && a3 > 0) {  // all positive, return min
            if (a1 > a2 || a1 > a3) {
                flag = true;  // changed if a1 is not min
            }
            double tmp = (a1 < a2) ? a1 : a2;
            return (tmp < a3) ? tmp : a3;
        }

        if (a1 < 0 && a2 < 0 && a3 < 0) {  // all negative, return max
            if (a1 < a2 || a1 < a3) {
                flag = true;  // changed if a1 is not max
            }
            double tmp = (a1 > a2) ? a1 : a2;
            return (tmp > a3) ? tmp : a3;
        }

        // return 0
        flag = true;  // changed
        return 0;
    }

    static void DG_recover(std::vector<double> &ret_u, size_t id_start,
                           size_t DG_k, double mean, double left,
                           double right) {
        if (DG_k == 0) { ret_u[id_start] = mean; }
        else if (DG_k == 1) {
            ret_u[id_start] = mean;
            ret_u[id_start + 1] = mean - left;
        }
        else if (DG_k >= 2) {
            ret_u[id_start] = mean;
            ret_u[id_start + 1] = -0.5 * left + 0.5 * right;
            ret_u[id_start + 2] = -mean + 0.5 * left + 0.5 * right;

            for (size_t j = 3; j <= DG_k; j++) {
                ret_u[id_start + j] = 0;  // ignore
            }
        }

        return;
    }

private:
    double m_tvb_M{0};
};
}  // namespace flux

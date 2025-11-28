#pragma once

#include <limits>
#include <string>

#include "expected.hpp"
#include "requires.h"

namespace flux::solver_deducing {

template <VarRequirements VarType, typename ExType>
class Solver {
public:
    auto run(this const auto &self, VarType var, ExType &ex, double t0,
             double tend) -> flux::expected<VarType, std::string> {
        if (tend <= t0) return var;

        double t = 0;
        bool stop_flag = false;
        constexpr auto iter_max = std::numeric_limits<std::size_t>::max();
        for (size_t iter = 0; iter < iter_max && (!stop_flag); ++iter) {
            var = self.update(var, ex, t, stop_flag, tend);
        }
        if (!stop_flag) { return flux::unexpected{std::string{"Iteration exceeds"}}; }

        return var;
    }
};

template <VarRequirements VarType, typename ExType>
class EulerSolver : public Solver<VarType, ExType> {
public:
    VarType update(this const auto &self, const VarType &var, ExType &ex,
                   double &t, bool &stop_flag, double tend) {
        double dt = self.get_dt(var, ex, t);
        if (t + dt >= tend && t < tend) {
            dt = tend - t;
            stop_flag = true;
        }

        auto var_n = self.pre_process(var, ex, t);

        VarType var2 = var_n + dt * self.op_L(var_n, ex, t);

        var2 = self.post_process(var2, ex, t);

        t += dt;
        return var2;
    }

    VarType post_process(const VarType &var, ExType &ex, double t) const {
        return var;
    }

    VarType pre_process(const VarType &var, ExType &ex, double t) const {
        return var;
    }
};

template <VarRequirements VarType, typename ExType>
class RK3Solver : public Solver<VarType, ExType> {
public:
    VarType update(this const auto &self, const VarType &var, ExType &ex,
                   double &t, bool &stop_flag, double tend) {
        double dt = self.get_dt(var, ex, t);
        if (t + dt >= tend && t < tend) {
            dt = tend - t;
            stop_flag = true;
        }

        auto var_n = self.pre_process(var, ex, t);

        VarType var1 = var_n + dt * self.op_L(var_n, ex, t);

        var1 = self.post_process_rk_stage(var1, ex, t);

        VarType var2 = (3.0 / 4) * var_n
                       + (1.0 / 4) * (var1 + dt * self.op_L(var1, ex, t + dt));

        var2 = self.post_process_rk_stage(var2, ex, t + dt);

        VarType var3 =
            (1.0 / 3) * var_n
            + (2.0 / 3) * (var2 + dt * self.op_L(var2, ex, t + dt / 2));

        var3 = self.post_process_rk_stage(var3, ex, t + dt / 2);

        var3 = self.post_process(var3, ex, t + dt);

        t += dt;
        return var3;
    }

    VarType post_process(const VarType &var, ExType &ex, double t) const {
        return var;
    }

    VarType pre_process(const VarType &var, ExType &ex, double t) const {
        return var;
    }

    VarType post_process_rk_stage(const VarType &var, ExType &ex,
                                  double t) const {
        return var;
    }
};
}  // namespace flux::solver_deducing

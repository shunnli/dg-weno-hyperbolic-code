#pragma once

#include <limits>
#include <string>

#include "expected.hpp"
#include "requires.h"

namespace flux::solver_crtp {

template <VarRequirements VarType, typename ExType, typename Derived>
class Solver {
public:
    auto run(VarType var, ExType &ex, double t0,
             double tend) const -> flux::expected<VarType, std::string> {
        if (tend <= t0) return var;

        double t = 0;
        bool stop_flag = false;
        constexpr auto iter_max = std::numeric_limits<std::size_t>::max();
        for (size_t iter = 0; iter < iter_max && (!stop_flag); ++iter) {
            var = derived().update(var, ex, t, stop_flag, tend);
        }
        if (!stop_flag) { return flux::unexpected{std::string{"Iteration exceeds"}}; }

        return var;
    }

protected:
    constexpr const Derived &derived() const {
        return static_cast<const Derived &>(*this);
    }
};

template <VarRequirements VarType, typename ExType, typename Derived>
class EulerSolver : public Solver<VarType, ExType, Derived> {
public:
    VarType update(const VarType &var, ExType &ex, double &t, bool &stop_flag,
                   double tend) const {
        double dt = derived().get_dt(var, ex, t);
        if (t + dt >= tend && t < tend) {
            dt = tend - t;
            stop_flag = true;
        }

        auto var_n = derived().pre_process(var, ex, t);

        VarType var2 = var_n + dt * derived().op_L(var_n, ex, t);

        var2 = derived().post_process(var2, ex, t);

        t += dt;
        return var2;
    }

    VarType post_process(const VarType &var, ExType &ex, double t) const {
        return var;
    }

    VarType pre_process(const VarType &var, ExType &ex, double t) const {
        return var;
    }

protected:
    constexpr const Derived &derived() const {
        return static_cast<const Derived &>(*this);
    }
};

template <VarRequirements VarType, typename ExType, typename Derived>
class RK3Solver : public Solver<VarType, ExType, Derived> {
public:
    VarType update(const VarType &var, ExType &ex, double &t, bool &stop_flag,
                   double tend) const {
        double dt = derived().get_dt(var, ex, t);
        if (t + dt >= tend && t < tend) {
            dt = tend - t;
            stop_flag = true;
        }

        auto var_n = derived().pre_process(var, ex, t);

        VarType var1 = var_n + dt * derived().op_L(var_n, ex, t);

        var1 = derived().post_process_rk_stage(var1, ex, t);

        VarType var2 =
            (3.0 / 4) * var_n
            + (1.0 / 4) * (var1 + dt * derived().op_L(var1, ex, t + dt));

        var2 = derived().post_process_rk_stage(var2, ex, t + dt);

        VarType var3 =
            (1.0 / 3) * var_n
            + (2.0 / 3) * (var2 + dt * derived().op_L(var2, ex, t + dt / 2));

        var3 = derived().post_process_rk_stage(var3, ex, t + dt / 2);

        var3 = derived().post_process(var3, ex, t + dt);

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

protected:
    constexpr const Derived &derived() const {
        return static_cast<const Derived &>(*this);
    }
};

}  // namespace flux::solver_crtp

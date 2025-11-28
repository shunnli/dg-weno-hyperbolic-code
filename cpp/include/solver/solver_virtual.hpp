#pragma once

#include <limits>
#include <string>

#include "expected.hpp"
#include "requires.h"

namespace flux::solver_virtual {

template <VarRequirements VarType, typename ExType>
class Solver {
public:
    auto run(VarType var, ExType &ex, double t0,
             double tend) const -> flux::expected<VarType, std::string> {
        if (tend <= t0) return var;

        double t = 0;
        bool stop_flag = false;
        constexpr auto iter_max = std::numeric_limits<std::size_t>::max();
        for (size_t iter = 0; iter < iter_max && (!stop_flag); ++iter) {
            var = update(var, ex, t, stop_flag, tend);
        }
        if (!stop_flag) {
            return flux::unexpected{std::string{"Iteration exceeds"}};
        }

        return var;
    }

    virtual VarType update(const VarType &var, ExType &ex, double &t,
                           bool &stop_flag, double tend) const = 0;

    virtual ~Solver() = default;
};

template <VarRequirements VarType, typename ExType>
class EulerSolver : public Solver<VarType, ExType> {
public:
    virtual double get_dt(const VarType &var, ExType &ex, double t) const = 0;

    virtual VarType op_L(const VarType &var, ExType &ex, double t) const = 0;

    virtual VarType post_process(const VarType &var, ExType &ex,
                                 double t) const {
        return var;
    }

    virtual VarType pre_process(const VarType &var, ExType &ex,
                                double t) const {
        return var;
    }

    VarType update(const VarType &var, ExType &ex, double &t, bool &stop_flag,
                   double tend) const override {
        double dt = get_dt(var, ex, t);
        if (t + dt >= tend && t < tend) {
            dt = tend - t;
            stop_flag = true;
        }

        auto var_n = this->pre_process(var, ex, t);

        VarType var2 = var_n + dt * op_L(var_n, ex, t);

        var2 = this->post_process(var2, ex, t);

        t += dt;
        return var2;
    }
};

template <VarRequirements VarType, typename ExType>
class RK3Solver : public Solver<VarType, ExType> {
public:
    virtual double get_dt(const VarType &var, ExType &ex, double t) const = 0;

    virtual VarType op_L(const VarType &var, ExType &ex, double t) const = 0;

    virtual VarType post_process(const VarType &var, ExType &ex,
                                 double t) const {
        return var;
    }

    virtual VarType pre_process(const VarType &var, ExType &ex,
                                double t) const {
        return var;
    }

    virtual VarType post_process_rk_stage(const VarType &var, ExType &ex,
                                          double t) const {
        return var;
    }

    VarType update(const VarType &var, ExType &ex, double &t, bool &stop_flag,
                   double tend) const override {
        double dt = get_dt(var, ex, t);
        if (t + dt >= tend && t < tend) {
            dt = tend - t;
            stop_flag = true;
        }

        auto var_n = this->pre_process(var, ex, t);

        VarType var1 = var_n + dt * op_L(var_n, ex, t);

        var1 = post_process_rk_stage(var1, ex, t);

        VarType var2 = (3.0 / 4) * var_n
                       + (1.0 / 4) * (var1 + dt * op_L(var1, ex, t + dt));

        var2 = post_process_rk_stage(var2, ex, t + dt);

        VarType var3 = (1.0 / 3) * var_n
                       + (2.0 / 3) * (var2 + dt * op_L(var2, ex, t + dt / 2));

        var3 = post_process_rk_stage(var3, ex, t + dt / 2);

        var3 = this->post_process(var3, ex, t + dt);

        t += dt;
        return var3;
    }
};
}  // namespace flux::solver_virtual

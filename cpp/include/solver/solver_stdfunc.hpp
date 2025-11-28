#pragma once

#include <functional>
#include <limits>
#include <string>

#include "expected.hpp"
#include "requires.h"

namespace flux::solver_stdfunc {

template <VarRequirements VarType, typename ExType>
class Solver {
public:
    using UpdateFunc = std::function<VarType(const VarType &, ExType &,
                                             double &, bool &, double)>;

    Solver &set_update(UpdateFunc update) {
        m_update = update;
        return *this;
    }

    auto run(VarType var, ExType &ex, double t0,
             double tend) const -> flux::expected<VarType, std::string> {
        if (m_update == nullptr) {
            return flux::unexpected{std::string{"update function is not set"}};
        }

        if (tend <= t0) return var;

        double t = 0;
        bool stop_flag = false;
        constexpr auto iter_max = std::numeric_limits<std::size_t>::max();
        for (size_t iter = 0; iter < iter_max && (!stop_flag); ++iter) {
            var = m_update(var, ex, t, stop_flag, tend);
        }
        if (!stop_flag) { return flux::unexpected{std::string{"Iteration exceeds"}}; }

        return var;
    }

protected:
    UpdateFunc m_update;
};

template <VarRequirements VarType, typename ExType>
class UpdaterFactory {
public:
    using OpFunc = std::function<VarType(const VarType &, ExType &, double)>;
    using DtFunc = std::function<double(const VarType &, ExType &, double)>;

    static auto get_euler_updater(OpFunc op_L, DtFunc get_dt,
                                  OpFunc pre_process, OpFunc post_process)
        -> Solver<VarType, ExType>::UpdateFunc {
        auto no_op = [](const VarType &var, ExType &ex, double t) {
            return var;
        };

        if (pre_process == nullptr) { pre_process = no_op; }
        if (post_process == nullptr) { post_process = no_op; }

        return [=](const VarType &var, ExType &ex, double &t, bool &stop_flag,
                   double tend) {
            double dt = get_dt(var, ex, t);
            if (t + dt >= tend && t < tend) {
                dt = tend - t;
                stop_flag = true;
            }
            auto var_n = pre_process(var, ex, t);

            VarType var2 = var_n + dt * op_L(var_n, ex, t);

            var2 = post_process(var2, ex, t);

            t += dt;
            return var2;
        };
    }

    static auto get_rk3_updater(
        OpFunc op_L, DtFunc get_dt, OpFunc pre_process, OpFunc post_process,
        OpFunc post_process_rk_stage) -> Solver<VarType, ExType>::UpdateFunc {
        auto no_op = [](const VarType &var, ExType &ex, double t) {
            return var;
        };

        if (pre_process == nullptr) { pre_process = no_op; }
        if (post_process == nullptr) { post_process = no_op; }
        if (post_process_rk_stage == nullptr) { post_process_rk_stage = no_op; }

        return [=](const VarType &var, ExType &ex, double &t, bool &stop_flag,
                   double tend) {
            double dt = get_dt(var, ex, t);
            if (t + dt >= tend && t < tend) {
                dt = tend - t;
                stop_flag = true;
            }

            auto var_n = pre_process(var, ex, t);

            VarType var1 = var_n + dt * op_L(var_n, ex, t);

            var1 = post_process_rk_stage(var1, ex, t);

            VarType var2 = (3.0 / 4) * var_n
                           + (1.0 / 4) * (var1 + dt * op_L(var1, ex, t + dt));

            var2 = post_process_rk_stage(var2, ex, t + dt);

            VarType var3 =
                (1.0 / 3) * var_n
                + (2.0 / 3) * (var2 + dt * op_L(var2, ex, t + dt / 2));

            var3 = post_process_rk_stage(var3, ex, t + dt / 2);

            var3 = post_process(var3, ex, t + dt);

            t += dt;
            return var3;
        };
    }
};
}  // namespace flux::solver_stdfunc

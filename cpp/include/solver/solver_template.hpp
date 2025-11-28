#pragma once

#include <limits>
#include <string>

#include "expected.hpp"
#include "requires.h"

namespace flux::solver_template {

template <typename UpdaterType, typename VarType, typename ExType>
concept UpdaterRequirements =
    requires(const UpdaterType &updater, const VarType &var, ExType &ex,
             double &t, bool &stop_flag, double tend) {
        { updater(var, ex, t, stop_flag, tend) } -> std::same_as<VarType>;
    };

template <VarRequirements VarType, typename ExType, typename UpdaterType>
    requires UpdaterRequirements<UpdaterType, VarType, ExType>
class Solver {
public:
    UpdaterType updater;

    auto run(VarType var, ExType &ex, double t0,
             double tend) const -> flux::expected<VarType, std::string> {
        if (tend <= t0) return var;

        double t = 0;
        bool stop_flag = false;
        constexpr auto iter_max = std::numeric_limits<std::size_t>::max();
        for (size_t iter = 0; iter < iter_max && (!stop_flag); ++iter) {
            var = updater(var, ex, t, stop_flag, tend);
        }
        if (!stop_flag) { return flux::unexpected{std::string{"Iteration exceeds"}}; }

        return var;
    }
};

template <typename OpType, typename VarType, typename ExType>
concept OpRequirements =
    requires(const OpType &op, const VarType &var, ExType &ex, double t) {
        { op(var, ex, t) } -> std::same_as<VarType>;
    };

template <typename GetDtType, typename VarType, typename ExType>
concept GetDtRequirements =
    requires(GetDtType get_dt, const VarType &var, ExType &ex, double t) {
        { get_dt(var, ex, t) } -> std::same_as<double>;
    };

template <VarRequirements VarType, typename ExType>
struct OpNull {
    auto operator()(const VarType &var, ExType &ex, double t) const {
        return var;
    }
};

template <VarRequirements VarType, typename ExType, typename OpType,
          typename GetDtType, typename PreProcessType, typename PostProcessType>
    requires OpRequirements<OpType, VarType, ExType>
             && GetDtRequirements<GetDtType, VarType, ExType>
             && OpRequirements<PreProcessType, VarType, ExType>
             && OpRequirements<PostProcessType, VarType, ExType>
class EulerUpdater {
public:
    OpType op_L;
    GetDtType get_dt;
    PreProcessType pre_process;
    PostProcessType post_process;

    VarType operator()(const VarType &var, ExType &ex, double &t,
                       bool &stop_flag, double tend) const {
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
    }
};

template <VarRequirements VarType, typename ExType, typename OpType,
          typename GetDtType, typename PreProcessType, typename PostProcessType,
          typename PostProcessRKStageType>
    requires OpRequirements<OpType, VarType, ExType>
             && GetDtRequirements<GetDtType, VarType, ExType>
             && OpRequirements<PreProcessType, VarType, ExType>
             && OpRequirements<PostProcessType, VarType, ExType>
             && OpRequirements<PostProcessRKStageType, VarType, ExType>
class RK3Updater {
public:
    OpType op_L;
    GetDtType get_dt;
    PreProcessType pre_process;
    PostProcessType post_process;
    PostProcessRKStageType post_process_rk_stage;

    VarType operator()(const VarType &var, ExType &ex, double &t,
                       bool &stop_flag, double tend) const {
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

        VarType var3 = (1.0 / 3) * var_n
                       + (2.0 / 3) * (var2 + dt * op_L(var2, ex, t + dt / 2));

        var3 = post_process_rk_stage(var3, ex, t + dt / 2);

        var3 = post_process(var3, ex, t + dt);

        t += dt;
        return var3;
    }
};
}  // namespace flux::solver_template

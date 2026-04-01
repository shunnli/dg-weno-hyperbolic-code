# FV-WENO5 and FD-WENO5 for 1D Burgers Equation

using LinearAlgebra
using Printf
using Plots

# ============================================================
# Gauss-Legendre quadrature (Golub-Welsch algorithm)
# ============================================================
function gauss_legendre(n)
    beta = [i / sqrt(4.0 * i^2 - 1) for i in 1:n-1]
    J = diagm(1 => beta, -1 => beta)
    F = eigen(Symmetric(J))
    perm = sortperm(F.values)
    nodes = F.values[perm]
    weights = 2.0 .* F.vectors[1, perm] .^ 2
    return nodes, weights
end

# ============================================================
# Burgers equation flux and flux splitting
# ============================================================
burgers_f(u) = u .^ 2 ./ 2
burgers_df(u) = u

function burgers_fhat_godunov(ul, ur)
    return ifelse.(
        ul .<= ur,
        ifelse.(ul .* ur .> 0, min.(burgers_f(ul), burgers_f(ur)), 0.0),
        max.(burgers_f(ul), burgers_f(ur)),
    )
end

function get_LF_flux(f, get_alpha)
    return (ul, ur) -> (f(ul) .+ f(ur)) ./ 2 .- get_alpha(ul, ur) ./ 2 .* (ur .- ul)
end

burgers_get_alpha(ul, ur) = max.(abs.(ul), abs.(ur))

burgers_f_split_l(u, c) = u .^ 2 ./ 4 .+ c ./ 2 .* u
burgers_f_split_r(u, c) = u .^ 2 ./ 4 .- c ./ 2 .* u

# ============================================================
# Exact Burgers solver: sinusoidal initial condition
# u(x,0) = a + b*sin(w*x + phi)
# ============================================================
function burgers_sin_solver_kernel(x, t)
    u = x ./ (π / 2 + t)
    for _ in 1:100_000
        du = (u .- sin.(x .- u .* t)) ./ (1 .+ cos.(x .- u .* t) .* t)
        u = u .- du
        maximum(abs.(du)) < 1e-10 && break
    end
    return u
end

function burgers_sin_solver(x, t, a, b, w, phi)
    x_tilde = w .* x .+ phi .- a .* w .* t
    x_tilde = mod.(x_tilde .+ π, 2π) .- π    # map to [-π, π]
    t_tilde = b * w * t
    return a .+ b .* burgers_sin_solver_kernel(x_tilde, t_tilde)
end

# Cell-averaged exact solution using 5-point Gauss quadrature
function burgers_sin_solver_mean(x, t, dx, a, b, w, phi)
    nodes, weights = gauss_legendre(5)
    x_nodes = x .+ (dx / 2) .* nodes'          # (n, 5) via broadcasting
    u_vals = burgers_sin_solver(x_nodes, t, a, b, w, phi)  # (n, 5)
    return 0.5 .* (u_vals * weights)            # (n,)
end

burgers_init_mean(x, dx, a, b) = burgers_sin_solver_mean(x, 0.0, dx, a, b, 1.0, 0.0)
burgers_init(x, a, b) = a .+ b .* sin.(x)

# ============================================================
# WENO5 reconstruction
# Returns (ul_hat, ur_hat):
#   ul_hat[i] = u_{i-1/2}^+  (left boundary of cell i, from inside)
#   ur_hat[i] = u_{i+1/2}^-  (right boundary of cell i, from inside)
# ============================================================
function compute_weno5_flux(u; do_left=true, do_right=true)
    EPSILON = 1e-6
    d_left = [3 / 10, 3 / 5, 1 / 10]
    d_right = [1 / 10, 3 / 5, 3 / 10]

    ul2 = circshift(u, 2)
    ul1 = circshift(u, 1)
    ur1 = circshift(u, -1)
    ur2 = circshift(u, -2)

    beta1 = @. (13 / 12) * (ul2 - 2 * ul1 + u)^2 + (1 / 4) * (ul2 - 4 * ul1 + 3 * u)^2
    beta2 = @. (13 / 12) * (ul1 - 2 * u + ur1)^2 + (1 / 4) * (ul1 - ur1)^2
    beta3 = @. (13 / 12) * (u - 2 * ur1 + ur2)^2 + (1 / 4) * (3 * u - 4 * ur1 + ur2)^2
    beta = hcat(beta1, beta2, beta3)

    ul_hat = zeros(length(u))
    ur_hat = zeros(length(u))

    if do_left
        p1 = @. -(1 / 6) * ul2 + (5 / 6) * ul1 + (1 / 3) * u
        p2 = @. (1 / 3) * ul1 + (5 / 6) * u - (1 / 6) * ur1
        p3 = @. (11 / 6) * u - (7 / 6) * ur1 + (1 / 3) * ur2
        ul_pre = hcat(p1, p2, p3)

        w_left = d_left' ./ (beta .+ EPSILON) .^ 2
        w_left = w_left ./ sum(w_left, dims=2)
        ul_hat = vec(sum(w_left .* ul_pre, dims=2))
    end

    if do_right
        p1 = @. (1 / 3) * ul2 - (7 / 6) * ul1 + (11 / 6) * u
        p2 = @. -(1 / 6) * ul1 + (5 / 6) * u + (1 / 3) * ur1
        p3 = @. (1 / 3) * u + (5 / 6) * ur1 - (1 / 6) * ur2
        ur_pre = hcat(p1, p2, p3)

        w_right = d_right' ./ (beta .+ EPSILON) .^ 2
        w_right = w_right ./ sum(w_right, dims=2)
        ur_hat = vec(sum(w_right .* ur_pre, dims=2))
    end

    return ul_hat, ur_hat
end

# ============================================================
# SSP-RK3 time integrator
# ============================================================
function SSP_RK3(u, dt, op_L)
    u1 = u .+ dt .* op_L(u)
    u2 = 3 / 4 .* u .+ 1 / 4 .* (u1 .+ dt .* op_L(u1))
    return 1 / 3 .* u .+ 2 / 3 .* (u2 .+ dt .* op_L(u2))
end

# ============================================================
# FV-WENO5 solver
# ============================================================
function FVWENO5(uh0, dx, tend, fhat, df)
    function op_L(u)
        ul_hat, ur_hat = compute_weno5_flux(u)
        fhat_r = fhat(ur_hat, circshift(ul_hat, -1))
        fhat_l = fhat(circshift(ur_hat, 1), ul_hat)
        return -1 / dx .* (fhat_r .- fhat_l)
    end

    uh = copy(uh0)
    tnow = 0.0
    while tnow < tend
        c = maximum(abs.(df(uh)))
        dtnow = min(1 / (2c) * dx^(5 / 3), tend - tnow)
        uh = SSP_RK3(uh, dtnow, op_L)
        tnow += dtnow
    end
    return uh
end

# ============================================================
# FD-WENO5 solver (flux-splitting form)
# ============================================================
function FDWENO5(uh0, dx, tend, fl, fr, df)
    function op_L(u)
        c = maximum(abs.(df(u)))
        hl_list = fl(u, c)
        hr_list = fr(u, c)
        _, hlr_hat = compute_weno5_flux(hl_list, do_left=false)
        hrl_hat, _ = compute_weno5_flux(hr_list, do_right=false)
        v = (hlr_hat .+ circshift(hrl_hat, -1)) .- (circshift(hlr_hat, 1) .+ hrl_hat)
        return -1 / dx .* v
    end

    uh = copy(uh0)
    tnow = 0.0
    while tnow < tend
        c = maximum(abs.(df(uh)))
        dtnow = min(1 / (2c) * dx^(5 / 3), tend - tnow)
        uh = SSP_RK3(uh, dtnow, op_L)
        tnow += dtnow
    end
    return uh
end

# ============================================================
# Convergence table utilities
# ============================================================
function convergence_order(xlist, ylist)
    return [-log(ylist[i+1] / ylist[i]) / log(xlist[i+1] / xlist[i]) for i in 1:length(xlist)-1]
end

function print_err_table(nums, err_l1, err_l2, err_inf)
    ord_l1 = convergence_order(nums, err_l1)
    ord_l2 = convergence_order(nums, err_l2)
    ord_linf = convergence_order(nums, err_inf)

    @printf("  %6s  %12s  %8s  %12s  %8s  %12s  %8s\n",
        "m", "L1 error", "L1 ord", "L2 error", "L2 ord", "Linf error", "Linf ord")
    for i in eachindex(nums)
        s1 = i == 1 ? "       -" : @sprintf("%8.4f", ord_l1[i-1])
        s2 = i == 1 ? "       -" : @sprintf("%8.4f", ord_l2[i-1])
        sinf = i == 1 ? "       -" : @sprintf("%8.4f", ord_linf[i-1])
        @printf("  %6d  %12.2e  %8s  %12.2e  %8s  %12.2e  %8s\n",
            nums[i], err_l1[i], s1, err_l2[i], s2, err_inf[i], sinf)
    end
end

# ============================================================
# Tests
# ============================================================
function FVWENO5_test()
    println("\n=== FV-WENO5: Convergence Test (t = 0.5) ===")
    tend = 0.5
    nums = [10, 20, 40, 80, 160, 320, 640]
    err_inf = zeros(length(nums))
    err_l1 = zeros(length(nums))
    err_l2 = zeros(length(nums))

    fhat = get_LF_flux(burgers_f, burgers_get_alpha)

    for (cnt, num) in enumerate(nums)
        dx = 2π / num
        x = collect(range(-π + dx / 2, step=dx, length=num))

        u0 = burgers_init_mean(x, dx, 0.5, 1.0)
        uh = FVWENO5(u0, dx, tend, fhat, burgers_df)
        u_exact = burgers_sin_solver_mean(x, tend, dx, 0.5, 1.0, 1.0, 0.0)

        err_inf[cnt] = maximum(abs.(u_exact .- uh))
        err_l1[cnt] = sum(abs.(u_exact .- uh)) * dx
        err_l2[cnt] = sqrt(sum((u_exact .- uh) .^ 2) * dx)
    end

    print_err_table(nums, err_l1, err_l2, err_inf)

    println("\n=== FV-WENO5: Post-shock solution (t = 1.5) ===")
    tend = 1.5
    num_ref = 320
    dx_ref = 2π / num_ref
    x_ref = collect(range(-π + dx_ref / 2, step=dx_ref, length=num_ref))
    u_ref = burgers_sin_solver_mean(x_ref, tend, dx_ref, 0.5, 1.0, 1.0, 0.0)

    for num in [40, 160]
        dx = 2π / num
        x = collect(range(-π + dx / 2, step=dx, length=num))
        u0 = burgers_init_mean(x, dx, 0.5, 1.0)
        uh = FVWENO5(u0, dx, tend, fhat, burgers_df)
        u_exact_on_coarse = burgers_sin_solver_mean(x, tend, dx, 0.5, 1.0, 1.0, 0.0)
        L1 = sum(abs.(u_exact_on_coarse .- uh)) * dx
        @printf("  n=%d: L1 error = %.2e\n", num, L1)
        p = plot(x_ref, u_ref; label="u_exact", lw=1, xlabel="x",
            title="FV-WENO5  t=$(tend)  n=$(num)", size=(600, 500), grid=true)
        scatter!(p, x, uh; label="u_h", markershape=:circle, ms=5, mc=:white, msw=1)
        savefig(p, joinpath(@__DIR__, "out/fv_weno5_shock_n$(num).png"))
    end
end

function FDWENO5_test()
    println("\n=== FD-WENO5: Convergence Test (t = 0.5) ===")
    tend = 0.5
    nums = [10, 20, 40, 80, 160, 320, 640]
    err_inf = zeros(length(nums))
    err_l1 = zeros(length(nums))
    err_l2 = zeros(length(nums))

    for (cnt, num) in enumerate(nums)
        dx = 2π / num
        x = collect(range(-π + dx / 2, step=dx, length=num))

        u0 = burgers_init(x, 0.5, 1.0)
        uh = FDWENO5(u0, dx, tend, burgers_f_split_l, burgers_f_split_r, burgers_df)
        u_exact = burgers_sin_solver(x, tend, 0.5, 1.0, 1.0, 0.0)

        err_inf[cnt] = maximum(abs.(u_exact .- uh))
        err_l1[cnt] = sum(abs.(u_exact .- uh)) * dx
        err_l2[cnt] = sqrt(sum((u_exact .- uh) .^ 2) * dx)
    end

    print_err_table(nums, err_l1, err_l2, err_inf)

    println("\n=== FD-WENO5: Post-shock solution (t = 1.5) ===")
    tend = 1.5
    num_ref = 320
    dx_ref = 2π / num_ref
    x_ref = collect(range(-π + dx_ref / 2, step=dx_ref, length=num_ref))
    u_ref = burgers_sin_solver(x_ref, tend, 0.5, 1.0, 1.0, 0.0)

    for num in [40, 160]
        dx = 2π / num
        x = collect(range(-π + dx / 2, step=dx, length=num))
        u0 = burgers_init(x, 0.5, 1.0)
        uh = FDWENO5(u0, dx, tend, burgers_f_split_l, burgers_f_split_r, burgers_df)
        u_exact_on_coarse = burgers_sin_solver(x, tend, 0.5, 1.0, 1.0, 0.0)
        L1 = sum(abs.(u_exact_on_coarse .- uh)) * dx
        @printf("  n=%d: L1 error = %.2e\n", num, L1)
        p = plot(x_ref, u_ref; label="u_exact", lw=1, xlabel="x",
            title="FD-WENO5  t=$(tend)  n=$(num)", size=(600, 500), grid=true)
        scatter!(p, x, uh; label="u_h", markershape=:circle, ms=5, mc=:white, msw=1)
        savefig(p, joinpath(@__DIR__, "out/fd_weno5_shock_n$(num).png"))
    end
end

# ============================================================
# Entry point
# ============================================================

mkpath(joinpath(@__DIR__, "out"))

FVWENO5_test()
FDWENO5_test()

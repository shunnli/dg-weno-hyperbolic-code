# Discontinuous Galerkin Method for 1D Scalar Conservation Laws
#
# Reference: Shu, Chi-Wang. "Discontinuous Galerkin Methods: General Approach
#            and Stability." (2008).

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
# Legendre basis functions on [-1, 1]
# ============================================================

# Evaluate P_0, ..., P_k at given points via three-term recurrence.
# Returns (k+1, n) matrix, or (k+1,) vector if points is a scalar.
function basis_eval(k, points)
    scalar_flag = isa(points, Real)
    pts = scalar_flag ? [Float64(points)] : collect(Float64.(reshape(points, :)))
    n = length(pts)

    result = zeros(k + 1, n)
    result[1, :] .= 1.0
    if k >= 1
        result[2, :] = pts
    end
    for deg in 1:k-1
        result[deg+2, :] = ((2deg + 1) .* pts .* result[deg+1, :] .- deg .* result[deg, :]) ./ (deg + 1)
    end
    return scalar_flag ? result[:, 1] : result
end

# Evaluate P'_0, ..., P'_k using the formula P'_n(x) = n/(x^2-1) * (x*P_n(x) - P_{n-1}(x)).
# Only valid at interior points (|x| < 1); used exclusively with Gauss quadrature nodes.
function basis_dx_eval(k, points)
    scalar_flag = isa(points, Real)
    pts = scalar_flag ? [Float64(points)] : collect(Float64.(reshape(points, :)))
    n = length(pts)

    result = zeros(k + 1, n)
    bvals = basis_eval(k, pts)     # (k+1, n)
    for deg in 1:k
        num = pts .* bvals[deg+1, :] .- bvals[deg, :]
        den = pts .^ 2 .- 1
        result[deg+1, :] = deg .* num ./ den
    end
    return scalar_flag ? result[:, 1] : result
end

# Analytical L2 inner products: integral_{-1}^{1} P_i(x)^2 dx = 2/(2i+1).
# Returns (k+1,) vector indexed from P_0 to P_k.
basis_inner_vec(k) = Float64[2 / (2 * (i - 1) + 1) for i in 1:k+1]

# ============================================================
# Burgers equation
# ============================================================
burgers_f(u) = u .^ 2 ./ 2
burgers_df(u) = u

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
    x_tilde = mod.(x_tilde .+ π, 2π) .- π     # map to [-π, π]
    t_tilde = b * w * t
    return a .+ b .* burgers_sin_solver_kernel(x_tilde, t_tilde)
end

burgers_init(x) = burgers_sin_solver(x, 0.0, 0.5, 1.0, 1.0, 0.0)
burgers_exact(x, t) = burgers_sin_solver(x, t, 0.5, 1.0, 1.0, 0.0)

# ============================================================
# Numerical flux: Lax-Friedrichs
# ============================================================
function get_LF_flux(f, get_alpha)
    return (ul, ur) -> (f(ul) .+ f(ur)) ./ 2 .- get_alpha(ul, ur) ./ 2 .* (ur .- ul)
end

burgers_get_alpha(ul, ur) = max.(abs.(ul), abs.(ur))

# ============================================================
# L2 projection onto piecewise Legendre polynomial space
# Returns DG coefficients (n, k+1) where uh[j, i] = c_{i-1}^{(j)}
# ============================================================
function l2_projection(f, x, dx, k, gk)
    nodes, weights = gauss_legendre(gk)
    M = basis_eval(k, nodes)       # (k+1, gk)
    inners = dx / 2 .* basis_inner_vec(k)  # (k+1,): denominator = dx/2 * 2/(2i+1)

    n = length(x)
    result = zeros(n, k + 1)
    for j in 1:n
        x_phys = x[j] .+ dx / 2 .* nodes   # (gk,) quadrature points in physical coords
        f_vals = f(x_phys)                   # (gk,)  — f must accept array input
        # Numerator: (dx/2) * sum_q f(x_q) * P_i(xi_q) * w_q  for each i
        result[j, :] = dx / 2 .* (M * (f_vals .* weights)) ./ inners
    end
    return result
end

# ============================================================
# Error computation: L∞, L1, L2 norms via Gauss quadrature
# ============================================================
function error4DG(uexact, uh, x, dx, t, k, gk)
    nodes, weights = gauss_legendre(gk)
    M = basis_eval(k, nodes)    # (k+1, gk); uh_val at quad pts = M' * uh[i,:]

    err_inf = 0.0
    err_l1 = 0.0
    err_l2 = 0.0
    for i in eachindex(x)
        u_ex = uexact(x[i] .+ dx / 2 .* nodes, t)   # (gk,)
        uh_val = M' * uh[i, :]                          # (gk,)
        diff = u_ex .- uh_val
        err_inf = max(err_inf, maximum(abs.(diff)))
        err_l1 += dot(abs.(diff), weights)
        err_l2 += dot(diff .^ 2, weights)
    end
    err_l1 *= dx / 2
    err_l2 = sqrt(dx / 2 * err_l2)
    return err_inf, err_l1, err_l2
end

# ============================================================
# TVD / TVB limiters
# ============================================================
function minmod(a1, a2, a3, _h)
    pos = (a1 .> 0) .& (a2 .> 0) .& (a3 .> 0)
    neg = (a1 .< 0) .& (a2 .< 0) .& (a3 .< 0)
    result = ifelse.(pos, min.(min.(a1, a2), a3), 0.0)
    result = ifelse.(neg, max.(max.(a1, a2), a3), result)
    return result
end

tvd() = minmod

function tvb(m)
    return (a1, a2, a3, h) -> ifelse.(abs.(a1) .< m .* h^2, a1, minmod(a1, a2, a3, h))
end

# ============================================================
# DG limiter post-processor
# Applies TVD/TVB limiting to DG coefficients while preserving
# cell averages. Degenerates to P2 for k >= 3.
# ============================================================
function dg_lim(lim, k, gk, h)
    @assert k >= 1

    nodes, weights = gauss_legendre(gk)
    M = basis_eval(k, nodes)        # (k+1, gk)
    vm = 0.5 .* (M * weights)        # (k+1,): projection of cell averages
    vl = basis_eval(k, -1.0)         # (k+1,): left  boundary basis values
    vr = basis_eval(k, 1.0)         # (k+1,): right boundary basis values

    # System matrix: rows = [cell-avg constraint, left-bdry, (right-bdry for k>=2)]
    if k == 1
        Mat = [vm[1] vm[2]
            vl[1] vl[2]]
    else
        Mat = [vm[1] vm[2] vm[3]
            vl[1] vl[2] vl[3]
            vr[1] vr[2] vr[3]]
    end

    function processer(uh)
        ul = uh * vl    # (n,): left  boundary values  u_{j-1/2}^+
        um = uh * vm    # (n,): cell averages           = uh[:, 1]
        ur = uh * vr    # (n,): right boundary values  u_{j+1/2}^-

        um_delta_r = circshift(um, -1) .- um
        um_delta_l = um .- circshift(um, 1)

        ul_mod = um .- lim(um .- ul, um_delta_l, um_delta_r, h)
        ur_mod = um .+ lim(ur .- um, um_delta_l, um_delta_r, h)

        n_cells, _ = size(uh)
        uh_mod = zeros(n_cells, k + 1)

        if k == 1
            B = [um'; ul_mod']               # (2, n)
            uh_mod = Matrix((Mat \ B)')      # (n, 2)
        else
            B = [um'; ul_mod'; ur_mod']      # (3, n)
            uh_mod[:, 1:3] = Matrix((Mat \ B)')  # (n, 3); higher coeffs remain 0
        end
        return uh_mod
    end

    return processer
end

# ============================================================
# SSP-RK3 with optional post-processing (e.g., slope limiter)
# ============================================================
function SSP_RK3_with_post_process(u, t, dt, op_L, post_process)
    u1 = post_process(u .+ dt .* op_L(u, t))
    u2 = post_process(3 / 4 .* u .+ 1 / 4 .* (u1 .+ dt .* op_L(u1, t + dt)))
    return post_process(1 / 3 .* u .+ 2 / 3 .* (u2 .+ dt .* op_L(u2, t + dt / 2)))
end

# ============================================================
# RKDG solver
#
# Arguments:
#   uh0   : initial DG coefficients (n, k+1)
#   dx    : cell width
#   tend  : final time
#   k     : polynomial degree
#   gk    : number of Gauss quadrature points (>= k+1)
#   f     : physical flux  f(u)
#   fhat  : numerical flux fhat(ul, ur)
#   df    : flux derivative f'(u) (for CFL estimation)
# Keyword:
#   lim   : limiter function (e.g. tvd(), tvb(5)); nothing = no limiter
# ============================================================
function RKDG(uh0, dx, tend, k, gk, f, fhat, df; lim=nothing)
    nodes, weights = gauss_legendre(gk)

    M = basis_eval(k, nodes)                         # (k+1, gk)
    D = (2 / dx) .* basis_dx_eval(k, nodes)         # (k+1, gk): d/dx in physical coords
    M_inv = (2 / dx) .* Diagonal(1.0 ./ basis_inner_vec(k))  # (k+1, k+1) inverse mass matrix

    vl = basis_eval(k, -1.0)   # (k+1,): basis at left  boundary
    vr = basis_eval(k, 1.0)   # (k+1,): basis at right boundary
    vc = basis_eval(k, 0.0)   # (k+1,): basis at cell center (for CFL estimate)

    function op_L(uh, ::Any)
        uh_l = circshift(uh, (1, 0))    # left  neighbor (shift rows down)
        uh_r = circshift(uh, (-1, 0))   # right neighbor (shift rows up)

        ul_m = uh_l * vr   # (n,): u at left  interface from left  cell
        ul_p = uh * vl   # (n,): u at left  interface from current cell
        ur_m = uh * vr   # (n,): u at right interface from current cell
        ur_p = uh_r * vl   # (n,): u at right interface from right cell

        fhat_l = fhat(ul_m, ul_p)   # (n,): numerical flux at left  interface
        fhat_r = fhat(ur_m, ur_p)   # (n,): numerical flux at right interface

        # Volume term: int_{I_j} f(u_h) * v_x dx
        tmp_u = uh * M                          # (n, gk): u at quadrature points
        tmp_fu = f.(tmp_u) .* weights'           # (n, gk): f(u) * weights (broadcast)
        int_fu_vx = (dx / 2) .* (tmp_fu * D')      # (n, k+1): volume integral

        # Boundary terms (outer products)
        fhatl_vl = fhat_l * vl'    # (n, k+1)
        fhatr_vr = fhat_r * vr'    # (n, k+1)

        # Residual: M_inv * (volume - right_flux + left_flux)
        delta = (int_fu_vx .- fhatr_vr .+ fhatl_vl) * M_inv
        return delta
    end

    post_process = (lim !== nothing && k > 0) ? dg_lim(lim, k, gk, dx) : identity

    # CFL-based base time step
    cfl = 1 / (2k + 1)
    dt0 = k > 2 ? cfl * dx^((k + 1) / 3) : cfl * dx

    uh = copy(uh0)
    t = 0.0
    while t < tend
        alpha_global = maximum(abs.(df(uh * vc)))
        dt = min(dt0 / alpha_global, tend - t)
        uh = SSP_RK3_with_post_process(uh, t, dt, op_L, post_process)
        t += dt
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
# Test functions
# ============================================================
function order_test(; init, exact, xleft, xright, nlist, tend, k, gk=7,
    f, fhat, df, lim=nothing)
    err_inf = zeros(length(nlist))
    err_l1 = zeros(length(nlist))
    err_l2 = zeros(length(nlist))

    for (cnt, n) in enumerate(nlist)
        dx = (xright - xleft) / n
        x = [xleft + (i - 0.5) * dx for i in 1:n]
        u0 = l2_projection(init, x, dx, k, gk)
        uh = RKDG(u0, dx, tend, k, gk, f, fhat, df; lim=lim)
        err_inf[cnt], err_l1[cnt], err_l2[cnt] = error4DG(exact, uh, x, dx, tend, k, gk)
    end
    print_err_table(nlist, err_l1, err_l2, err_inf)
end

function shock_test(; init, exact, xleft, xright, nlist, n_ref, tend, k, gk=7,
    f, fhat, df, lim=nothing, name::String="")
    scheme_name = isempty(name) ? "P$(k)" : name
    @printf("  (Reference: n_ref=%d exact solution at t=%.1f)\n", n_ref, tend)

    dx_ref = (xright - xleft) / n_ref
    x_ref = [xleft + (i - 0.5) * dx_ref for i in 1:n_ref]
    u_exact_ref = exact(x_ref, tend)

    vc = basis_eval(k, 0.0)    # evaluate DG polynomial at cell center

    for n in nlist
        dx = (xright - xleft) / n
        x = [xleft + (i - 0.5) * dx for i in 1:n]
        u0 = l2_projection(init, x, dx, k, gk)
        uh = RKDG(u0, dx, tend, k, gk, f, fhat, df; lim=lim)

        # Evaluate DG polynomial at cell centers
        uh_center = [dot(vc, uh[i, :]) for i in 1:n]
        u_ex = exact(x, tend)
        L1 = sum(abs.(u_ex .- uh_center)) * dx
        @printf("  n=%3d: L1 error = %.2e\n", n, L1)
        p = plot(x_ref, u_exact_ref; label="u_exact", lw=1, xlabel="x",
            title="DG $(scheme_name)  t=$(tend)  n=$(n)", size=(600, 500), grid=true)
        scatter!(p, x, uh_center; label="u_h", markershape=:circle, ms=3, mc=:white, msw=1)
        savefig(p, joinpath(@__DIR__, "out/dg_$(lowercase(scheme_name))_shock_n$(n).png"))
    end
end

# ============================================================
# Burgers equation characteristics / solution evolution plot
# (corresponds to Python notebook cell 5)
# ============================================================
function plot_burgers_evolution()
    dx = 2π / 100
    x = [dx / 2 + (i - 1) * dx for i in 1:100]
    p = plot(; xlabel="x", size=(600, 500), grid=true)
    for tend in [0.0, 0.5, 1.0, 2.0, 3.0]
        u = burgers_sin_solver(x, tend, 0.0, 1.0, 1.0, 0.0)
        if tend >= 1.0
            u0 = sin.(x)
            plot!(p, x .+ u0 .* tend, u0; ls=:dash, lw=1, label="")
        end
        plot!(p, x, u; lw=1, label="t=$(tend)")
    end
    savefig(p, joinpath(@__DIR__, "out/burgers_evolution.png"))
end

# ============================================================
# Main: run all tests
# ============================================================

mkpath(joinpath(@__DIR__, "out"))

lf_flux = get_LF_flux(burgers_f, burgers_get_alpha)

plot_burgers_evolution()

# --- DG P0 ---
println("=== DG P0: Convergence Test (t = 0.5) ===")
order_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[20, 40, 80, 160, 320], tend=0.5,
    k=0, f=burgers_f, df=burgers_df, fhat=lf_flux)
println("\n=== DG P0: Shock Test (t = 1.5) ===")
shock_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[40, 160], n_ref=320, tend=1.5,
    k=0, f=burgers_f, df=burgers_df, fhat=lf_flux)

# --- DG P1 ---
println("\n=== DG P1: Convergence Test (t = 0.5) ===")
order_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[20, 40, 80, 160, 320], tend=0.5,
    k=1, f=burgers_f, df=burgers_df, fhat=lf_flux)
println("\n=== DG P1: Shock Test (t = 1.5) ===")
shock_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[40, 160], n_ref=320, tend=1.5,
    k=1, f=burgers_f, df=burgers_df, fhat=lf_flux)

# --- DG P2 ---
println("\n=== DG P2: Convergence Test (t = 0.5) ===")
order_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[20, 40, 80, 160, 320], tend=0.5,
    k=2, f=burgers_f, df=burgers_df, fhat=lf_flux)
println("\n=== DG P2: Shock Test (t = 1.5) ===")
shock_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[40, 160], n_ref=320, tend=1.5,
    k=2, f=burgers_f, df=burgers_df, fhat=lf_flux)

# --- DG P2 + TVD ---
println("\n=== DG P2 + TVD: Convergence Test (t = 0.5) ===")
order_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[20, 40, 80, 160, 320], tend=0.5,
    k=2, f=burgers_f, df=burgers_df, fhat=lf_flux, lim=tvd())
println("\n=== DG P2 + TVD: Shock Test (t = 1.5) ===")
shock_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[40, 160], n_ref=320, tend=1.5,
    k=2, f=burgers_f, df=burgers_df, fhat=lf_flux, lim=tvd(), name="P2_TVD")

# --- DG P2 + TVB (M=5) ---
println("\n=== DG P2 + TVB (M=5): Convergence Test (t = 0.5) ===")
order_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[20, 40, 80, 160, 320], tend=0.5,
    k=2, f=burgers_f, df=burgers_df, fhat=lf_flux, lim=tvb(5))
println("\n=== DG P2 + TVB (M=5): Shock Test (t = 1.5) ===")
shock_test(init=burgers_init, exact=burgers_exact,
    xleft=-π, xright=π, nlist=[40, 160], n_ref=320, tend=1.5,
    k=2, f=burgers_f, df=burgers_df, fhat=lf_flux, lim=tvb(5), name="P2_TVB_M5")

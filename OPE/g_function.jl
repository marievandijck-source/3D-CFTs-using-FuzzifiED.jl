using JLD2, FuzzifiED, LinearAlgebra, Printf, Glob, SparseArrays
using Suppressor
include("Functions.jl")

# ─────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────

result_dir = "results/O3"
all_data   = load_all_results(result_dir)
sort!(all_data, by = d -> (d.nm, d.E))

# ─────────────────────────────────────────────
# Compute g_ε(h_c ± δh) for each system size
# ─────────────────────────────────────────────

results_table = []

for nm in [5, 6, 7, 8, 9, 10]
    println("nm = $nm")
    no  = 4 * nm
    nf  = 4
    h_c = 14.992 + 6.59 * nm ^ (-2.16 / 2)

    data_nm = filter(d -> d.nm == nm, all_data)

    # Ground state and reference eigenvectors
    eig_0  = data_nm[1].vec
    σ      = first(sort(select_states(0, 2, 1, data_nm), by = d -> d.Δ))
    ∂σ     = first(sort(select_states(2, 2, 1, data_nm), by = d -> d.Δ))

    # Hamiltonian parts for each symmetry sector
    @suppress begin
        global H_σ1,  H_σh  = H_parts(1, 1,  1, no, nm, nf)
        global H_∂σ1, H_∂σh = H_parts(1, 1, -1, no, nm, nf)
        global H_01,  H_0h  = H_parts(0, 1,  1, no, nm, nf)
    end

    for δh in [0.1, 0.05, 0.01, 0.005, 0.001]
        c_plus, g_plus = solve_c_gε(
            h_c + δh,
            H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh,
            eig_0, σ.vec, ∂σ.vec
        )
        c_min, g_min = solve_c_gε(
            h_c - δh,
            H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh,
            eig_0, σ.vec, ∂σ.vec
        )
        push!(results_table, (;
            nm, δh, h_c,
            g_plus, g_min,
            c_plus, c_min,
        ))
    end
end

# ─────────────────────────────────────────────
# Print results
# ─────────────────────────────────────────────

println("─"^75)
@printf "%-4s  %-8s  %-10s  %-12s  %-12s\n" "nm" "δh" "h_c" "g_plus" "g_min"
println("─"^75)

current_nm = -1
for r in results_table
    if r.nm != current_nm
        current_nm = r.nm
        println("─"^75)
    end
    @printf "%-4d  %-8.4f  %-10.4f  %-12.6f  %-12.6f\n" r.nm r.δh r.h_c r.g_plus r.g_min
end
println("─"^75)

# ─────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────

outdir  = "results"
isdir(outdir) || mkpath(outdir)
outfile = joinpath(outdir, "g_results.jld2")
@save outfile results_table
println("Saved to: ", outfile)

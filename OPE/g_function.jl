using JLD2
using Dates
using FuzzifiED
using LinearAlgebra
using Printf
using Glob
using SparseArrays



include("Functions.jl")

result_dir = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/HPC/O3"

all_data   = load_all_results(result_dir)

results_table = []

for nm = [5,6,7,8,9,10]
    println("start", nm)

    no = 4*nm
    nf = 4
    h_c = 14.992 + 6.59 * nm ^ (-2.16 / 2) 

    
    sort!(all_data, by = d -> (d.nm, d.E))
    data_nm = filter(d -> d.nm == nm, all_data)

    E_0 = data_nm[1].E

    state_σ = select_states(0,2,1,1,1, data_nm) 
    state_∂σ = select_states(2,2,1,1,-1, data_nm)

    σ  = first(sort(state_σ,  by = d -> d.Δ))
    ∂σ = first(sort(state_∂σ, by = d -> d.Δ))

    eig_σ_sec  = σ.vec
    eig_∂σ_sec = ∂σ.vec
    eig_0_sec = data_nm[1].vec

    H_σ1, H_σh = H_parts(1, 1, 1, no, nm, nf)
    H_∂σ1, H_∂σh = H_parts(1, 1, -1, no, nm, nf)
    H_01, H_0h = H_parts(0, 1, 1, no, nm, nf)
    

    for δh = [0.1, 0.05, 0.01, 0.005, 0.001]
        c_plus, g_plus = solve_c_gε_fast(h_c + δh, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0_sec, eig_σ_sec, eig_∂σ_sec)
        c_min, g_min = solve_c_gε_fast(h_c - δh, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0_sec, eig_σ_sec, eig_∂σ_sec)
        push!(results_table, (
            nm     = nm,
            δh     = δh,
            h_c    = h_c,
            g_plus = g_plus,
            g_min  = g_min,
            c_min = c_min,
            c_plus = c_plus,
        ))
    end 
end 

println("─"^75)
@printf "%-4s  %-8s  %-10s  %-12s  %-12s\n" "nm" "δh" "h_c" "g_plus" "g_min"
println("─"^75)

current_nm = -1
for r in results_table
    if r.nm != current_nm
        global current_nm = r.nm
        println("─"^75)
    end
    @printf "%-4d  %-8.4f  %-10.4f  %-12.6f  %-12.6f\n" r.nm r.δh r.h_c r.g_plus r.g_min
end
println("─"^75)

# ── Opslaan ───────────────────────────────────────────────────
outdir  = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/Julia/OPE/results"
isdir(outdir) || mkpath(outdir)
outfile = joinpath(outdir, "g_results.jld2")

@save outfile results_table
println("Opgeslagen in: ", outfile)
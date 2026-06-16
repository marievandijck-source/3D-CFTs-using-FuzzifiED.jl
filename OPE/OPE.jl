using JLD2, FuzzifiED, LinearAlgebra, Printf, Glob, SparseArrays
using Suppressor
include("Functions.jl")

# ─────────────────────────────────────────────
# Configuration
#
# Set the symmetry quantum numbers of the target state and the
# output filename. Available combinations are listed below.
#

file = "S1L0-.jld2"
l2   = [2]    # angular momentum: l(l+1) ∈ {0, 2, 6, 12, …}
c2   = [2]    # spin Casimir:     c(c+1) ∈ {0, 2, 6, 12, …}
I    = [1]    # index within symmetry sector
Z    = 1      # Z₂ parity (1 = odd, 0 = even)
X    = 1
R    = 1

# δh index: 1→0.1, 2→0.05, 3→0.01, 4→0.005, 5→0.001
δh_values = [0.1, 0.05, 0.01, 0.005, 0.001]
index     = 4
δh        = δh_values[index]

# ─────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────

result_dir = "results/O3"
all_data   = load_all_results(result_dir)
sort!(all_data, by = d -> (d.nm, d.E))

g_file = "results/g_results.jld2"
@load g_file results_table
sort!(results_table, by = d -> d.nm)

# ─────────────────────────────────────────────
# Define included system sizes
# ─────────────────────────────────────────────
sizes = [5, 6, 7, 8, 9, 10]

# ─────────────────────────────────────────────
# Compute OPE coefficients
# ─────────────────────────────────────────────

OPEs = []

for nm in sizes
    no  = 4 * nm
    nf  = 4
    h_c = 14.992 + 6.59 * nm ^ (-2.16 / 2)

    data_nm = filter(d -> d.nm == nm, all_data)
    g_data  = filter(d -> d.nm == nm, results_table)

    g_plus = g_data[index].g_plus
    g_min  = g_data[index].g_min
    c_plus = g_data[index].c_plus
    c_min  = g_data[index].c_min

    # Hamiltonian parts
    @suppress begin
        global H_01, H_0h = H_parts(0, 1, 1, no, nm, nf)
        global H_1,  H_h  = ([Z, X, R] == [0, 1, 1]) ?
                              (H_01, H_0h) :
                              H_parts(Z, X, R, no, nm, nf)
    end

    # Ground-state energy at h_c ± δh
    state_0 = data_nm[1].vec
    @suppress begin
        global E_0_p = real(state_0' * H_01 * state_0 + (h_c + δh) * state_0' * H_0h * state_0)
        global E_0_m = real(state_0' * H_01 * state_0 + (h_c - δh) * state_0' * H_0h * state_0)
    end

    for nr in eachindex(l2)
        states = select_states(l2[nr], c2[nr], Z, data_nm)
        state  = states[I[nr]]
        vec    = state.vec
        Δ      = state.Δ

        @suppress begin
            global E_p = real(vec' * H_1 * vec + (h_c + δh) * vec' * H_h * vec)
            global E_m = real(vec' * H_1 * vec + (h_c - δh) * vec' * H_h * vec)
        end

        f_ope = ((E_p - E_0_p) - (E_m - E_0_m) - (c_plus - c_min) * Δ) / (g_plus - g_min)

        println("nm=$nm  Δ=$(round(Δ; digits=6))  f=$(round(f_ope; digits=6))  " *
                "(Z=$Z, X=$X, R=$R, c2=$(c2[nr]), l2=$(l2[nr]))")

        push!(OPEs, (;
            nm, Z, X, R,
            l2    = l2[nr],
            c2    = c2[nr],
            index = I[nr],
            δh, h_c, Δ,
            f_ope,
        ))
    end
end

# ─────────────────────────────────────────────
# Print results
# ─────────────────────────────────────────────

println("─"^80)
@printf "%-4s  %-3s  %-3s  %-3s  %-4s  %-4s  %-8s  %-10s  %-12s\n" \
    "nm" "Z" "X" "R" "l2" "c2" "δh" "Δ" "f_ope"
println("─"^80)
for r in OPEs
    @printf "%-4d  %-3d  %-3d  %-3d  %-4d  %-4d  %-8.4f  %-10.6f  %-12.6f\n" \
        r.nm r.Z r.X r.R r.l2 r.c2 r.δh r.Δ r.f_ope
    println("─"^80)
end

# ─────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────

outdir  = "results"
isdir(outdir) || mkpath(outdir)
@save joinpath(outdir, file) OPEs
println("Saved to: ", joinpath(outdir, file))

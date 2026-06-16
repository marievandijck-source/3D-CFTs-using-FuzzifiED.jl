using FuzzifiED
using LinearAlgebra
using Glob
using JLD2

# ─────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────

"""
    load_all_results(result_dir)

Load all `O3_nm*.jld2` files from `result_dir` and return a flat list of named tuples,
one per eigenstate, with fields `(nm, Δ, E, l2, c2, Z, X, R, vec)`.
"""
function load_all_results(result_dir)
    files    = glob("O3_nm*.jld2", result_dir)
    all_data = []

    for f in sort(files)
        println("Loading: ", basename(f))
        try
            local nm, enrg_0, enrg_ϕ, results, eigenvectors
            @load f nm enrg_0 enrg_ϕ results eigenvectors

            for (key, res) in results
                Z, X, R = key
                vecs    = eigenvectors[key]   # columns = eigenvectors

                for i in eachindex(res.energy)
                    push!(all_data, (
                        nm  = nm,
                        Δ   = 0.518936 * (real(res.energy[i]) - enrg_0) / (enrg_ϕ - enrg_0),
                        E   = real(res.energy[i]),
                        l2  = real(res.l2[i]),
                        c2  = real(res.c2[i]),
                        Z   = Z,
                        X   = X,
                        R   = R,
                        vec = vecs[:, i],
                    ))
                end
            end
            println("  ✓ nm=$nm loaded")
        catch e
            @warn "Skipped: $(basename(f)) — $e"
        end
    end

    return all_data
end

# ─────────────────────────────────────────────
# State selection
# ─────────────────────────────────────────────

"""
    select_states(l2, c2, Z, data_nm; atol = 0.4)

Return all states in `data_nm` matching the given `l2`, `c2`, and `Z` quantum numbers.
"""
function select_states(l2, c2, Z, data_nm; atol = 0.4)
    return filter(d ->
        isapprox(d.l2, l2; atol) &&
        isapprox(d.c2, c2; atol) &&
        d.Z == Z,
        data_nm)
end

# ─────────────────────────────────────────────
# Hamiltonian construction
# ─────────────────────────────────────────────

"""
    H_parts(Z, X, R, no, nm, nf)

Construct and return the two operators `(H_1, H_h)` for the O(3) model, where
`H_1` is the interaction part and `H_h` is the perturbation (magnetic field) part.
"""
function H_parts(Z, X, R, no, nm, nf)
    mat_0 = diagm([0, 0, 0, 1.0])
    mat_V = [
        [ 0 0 0  1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0 -1 0 ] / √2,
        [ 0 0 0 -1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0  1 0 ] / √2 * im,
        [ 0 0 0  0 ; 0 0 0 1 ; 0 0 0  0 ; 0 1  0 0 ]
    ]
    qnd = [
        GetNeQNDiag(no),
        GetLz2QNDiag(nm, nf),
        GetFlavQNDiag(nm, nf, Dict([1 => 1, 3 => -1])),
        GetFlavQNDiag(nm, nf, Dict([f => 1 for f in 1:nf-1]), 0, 2)
    ]
    qnf = [
        GetRotyQNOffd(nm, nf),
        GetFlavPermQNOffd(nm, nf, Dict([1 => 3, 3 => 1]))
    ]

    cfs = Dict(Z_val => Confs(no, [nm, 0, 0, Z_val], qnd) for Z_val in 0:1)
    bs  = Basis(cfs[Z], [R, X], qnf)

    tms_1 = SimplifyTerms(
        GetDenIntTerms(nm, nf, [6.5, 1.0])
        - 1.4 * 2 * GetDenIntTerms(nm, nf, [6.5, 1.0], mat_V)
    )
    tms_h = SimplifyTerms(-2 * GetPolTerms(nm, nf, mat_0))

    return Operator(bs, tms_1), Operator(bs, tms_h)
end

# ─────────────────────────────────────────────
# g_ε and c/R solver
# ─────────────────────────────────────────────

"""
    solve_c_gε(h, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh,
               eig_0, eig_σ, eig_∂σ)

Solve a 2×2 linear system to extract `c/R` and `g_ε(h)` from the energies of the
ground state, the σ state, and the ∂σ state at coupling `h`.

Returns `(c_over_R, g_ε)`.
"""
function solve_c_gε(h, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0, eig_σ, eig_∂σ)
    Δ_σ   = 0.518936
    Δ_ε   = 1.59488
    f_σεσ = 0.525
    A_σε  = 1 + Δ_ε * (Δ_ε - 3) / (6 * Δ_σ)

    E_0  = real(eig_0'  * H_01  * eig_0  + h * eig_0'  * H_0h  * eig_0)
    E_σ  = real(eig_σ'  * H_σ1  * eig_σ  + h * eig_σ'  * H_σh  * eig_σ)  - E_0
    E_∂σ = real(eig_∂σ' * H_∂σ1 * eig_∂σ + h * eig_∂σ' * H_∂σh * eig_∂σ) - E_0

    M = [Δ_σ       f_σεσ      ;
         Δ_σ + 1   f_σεσ * A_σε]

    x = M \ [E_σ, E_∂σ]
    return x[1], x[2]   # c/R, g_ε
end


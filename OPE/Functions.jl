function load_all_results(result_dir)
    files = glob("O3_nm*.jld2", result_dir)
    all_data = []

    for f in sort(files)
        println("Laden: ", basename(f))

        try
            local nm, enrg_0, enrg_ϕ, results, eigenvectors
            @load f nm enrg_0 enrg_ϕ results eigenvectors

            for (key, res) in results
                Z, X, R = key
                vecs = eigenvectors[key]  # matrix: dim_sector × n_eigenvalues

                for i in eachindex(res.energy)
                    E   = real(res.energy[i])
                    l2  = real(res.l2[i])
                    c2  = real(res.c2[i])
                    Δ   = 0.518936 * (E - enrg_0) / (enrg_ϕ - enrg_0)
                    vec = vecs[:, i]  # kolom i = i-de eigenvector van deze sector

                    push!(all_data, (
                        nm  = nm,
                        Δ   = Δ,
                        E   = E,
                        l2  = l2,
                        c2  = c2,
                        Z   = Z,
                        X   = X,
                        R   = R,
                        vec = vec,
                    ))
                end
            end

            println("  ✓ nm=$nm geladen")

        catch e
            @warn "Overgeslagen: $(basename(f)) — $e"
            showerror(stdout, e)
            println()
        end
    end
    
    return all_data
end


function select_states(l2,c2,Z, data_nm)
    state = filter(d ->
    isapprox(d.l2, l2, atol=0.4) &&
    isapprox(d.c2, c2, atol=0.4) &&
    d.Z == Z,
    data_nm)
    return state
end

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
        GetFlavQNDiag(nm, nf, Dict([f => 1 for f = 1 : nf - 1]), 0, 2)
    ]

    qnf = [
        GetRotyQNOffd(nm, nf),
        GetFlavPermQNOffd(nm, nf, Dict([1 => 3, 3 => 1]))
    ]

    cfs = Dict{Int64, Confs}()
    for Z = 0 : 1
        cfs[Z] = Confs(no, [nm, 0, 0, Z], qnd)
    end

    bs = Basis(cfs[Z], [R, X], qnf)
    
    tms_0 = SimplifyTerms(
        GetDenIntTerms(nm, nf, [6.5, 1.0])
        - 1.4 * 2 * GetDenIntTerms(nm, nf, [6.5, 1.0], mat_V)
    )
    tms_h = SimplifyTerms(-2 * GetPolTerms(nm, nf, mat_0))

    op_0 = Operator(bs, tms_0)
    op_h = Operator(bs, tms_h)
    return op_0, op_h
end



function solve_c_gε_fast(h_c, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0_sec, eig_σ_sec, eig_∂σ_sec)
    Δ_σ  = 0.518936
    Δ_ε  = 1.59488 
    f_σεσ = 0.525
    A_σε  = 1 + Δ_ε * (Δ_ε - 3) / (6 * Δ_σ)


    E_0 = (eig_0_sec' * H_01 * eig_0_sec) + h_c * (eig_0_sec' * H_0h * eig_0_sec)
    println("groundstate= ", E_0)
    E_σ = (eig_σ_sec' * H_σ1 * eig_σ_sec) + h_c * (eig_σ_sec' * H_σh * eig_σ_sec) - E_0
    println("sigma= ",E_σ)

    println("sum", E_0 + E_0 + E_σ)
    E_∂σ = (eig_∂σ_sec' * H_∂σ1 * eig_∂σ_sec) + h_c * (eig_∂σ_sec'  * H_∂σh * eig_∂σ_sec) - E_0

    M = [Δ_σ      f_σεσ ;
         Δ_σ + 1  f_σεσ * A_σε]
    b = [E_σ, E_∂σ]
    x = M \ b
    println("x=",x[1])
    return x[1], x[2]  # c/R, g_ε
end

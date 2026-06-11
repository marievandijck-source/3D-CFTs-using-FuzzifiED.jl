using JLD2
using Printf
using Glob

function load_all_results(result_dir)
    files = glob("O3_nm*.jld2", result_dir)
    all_data = []

    for f in sort(files)
        println("Laden: ", basename(f))

        try
            local nm, enrg_0, enrg_ϕ, results
            @load f nm enrg_0 enrg_ϕ results

            for (key, res) in results
                Z, X, R = key

                for i in eachindex(res.energy)
                    E  = real(res.energy[i])
                    l2 = real(res.l2[i])
                    c2 = real(res.c2[i])
                    Δ  = 0.518936 * (E - enrg_0) / (enrg_ϕ - enrg_0)

                    push!(all_data, (
                        nm = nm,
                        Δ  = Δ,
                        E  = E,
                        l2 = l2,
                        c2 = c2,
                        Z  = Z,
                        X  = X,
                        R  = R,
                    ))
                end
            end

            println("  ✓ nm=$nm geladen")

        catch e
            @warn "Overgeslagen: $(basename(f)) — $e"
        end
    end

    return all_data
end

# ── Laad alle bestanden ───────────────────────────────────────
result_dir = joinpath(@__DIR__, "O3")
all_data = load_all_results(result_dir)

# ── Sorteer op nm en energie ──────────────────────────────────
sort!(all_data, by = d -> (d.nm, d.E))


function print_results(all_data)
    current_nm = -1
    for d in all_data
        if d.nm != current_nm
            current_nm = d.nm
            println("\n===== nm = $current_nm =====")
            @printf "%-8s  %-8s  %-8s  %-4s  %-4s  %-4s\n" "Δ" "L²" "C₂" "Z" "X" "R"
            println("─"^45)
        end
        @printf "%-8.4f  %-8.3f  %-8.3f  %-4d  %-4d  %-4d\n" d.Δ d.l2 d.c2 d.Z d.X d.R
    end
end

result_dir = joinpath(@__DIR__, "O3")
all_data = load_all_results(result_dir)
sort!(all_data, by = d -> (d.nm, d.E))



# ── Opslaan ───────────────────────────────────────────────────
outfile = joinpath(result_dir, "O3_combined2.jld2")
@save outfile all_data
println("\nOpgeslagen in $outfile")

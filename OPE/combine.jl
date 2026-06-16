using JLD2

# ─────────────────────────────────────────────
# Combine all OPE result files into one .jld2 file.
# Set `result_dir` to the folder containing the individual output files.
# ─────────────────────────────────────────────

result_dir = "results"

all_OPEs = []
for f in filter(readdir(result_dir; join = true)) do f
        endswith(f, ".jld2") && basename(f) != "combined_OPEs.jld2"
    end
    append!(all_OPEs, load(f)["OPEs"])
end

println("Total OPE entries: ", length(all_OPEs))

outfile = joinpath(result_dir, "combined_OPEs.jld2")
@save outfile all_OPEs
println("Saved to: ", outfile)

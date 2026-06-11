using JLD2

result_dir = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/julia/OPE/results/L2S1+"

all_OPEs = []
for f in filter(f -> endswith(f, ".jld2") && f != "combined_OPEs.jld2", readdir(result_dir, join=true))
    local d = load(f)
    append!(all_OPEs, d["OPEs"])
end

println("Totaal aantal OPEs: ", length(all_OPEs))

outfile = joinpath(result_dir, "combined_OPEs.jld2")
@save outfile all_OPEs
println("Opgeslagen in: ", outfile)
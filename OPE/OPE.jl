#Needed: Functions.jl
#Needed: Data for states and date for g

#L0S0- : ok
#L0S0+ : ok
#L1S0- : ok
#L1S0+ : ok
#L2S0- : ok
#L2S0+ : ok

#L0S1- :
#L0S1+
#L1S1- : ok
#L1S1+ : ok
#L2S1-
#L2S1+

#L0S2-
#L0S2+
#L1S2-
#L1S2+
#L2S2-
#L2S2+

using JLD2
using Dates
using FuzzifiED
using LinearAlgebra
using Printf
using Glob
using SparseArrays
using Suppressor

#Symmetry 
file = "S1L0n2.jld2"
l2 = [2] #Values 0, 2, 6, 12, 20 ... (= l(l+1)) #To select state
c2 = [2] #Values 0, 2, 6, 12, 20 ... (= c(c+1)) #To select state
I = [1] #Order within symmetry sector 
Z = 1 #To construct H
X= 1 #To construct H
R= 1 #To construct H

#Selection of δh [0.1, 0.05, 0.01, 0.005, 0.001]
index = 4 #δh = 0.001
δh = [0.1, 0.05, 0.01, 0.005, 0.001][index]

#Functions
include("Functions.jl")

#Data
result_dir = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/HPC/O3"
all_data   = load_all_results(result_dir)
sort!(all_data, by = d -> (d.nm, d.E))

OPEs = []

for nm in [5,6,7,8,9,10]

    no = 4*nm
    nf = 4
    h_c = 14.992 + 6.59 * nm ^ (-2.16 / 2)

    data_nm = filter(d -> d.nm == nm, all_data)

    #Data from the g
    g_file = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/julia/OPE/results/g_results.jld2"
    @load g_file results_table
    sort!(results_table, by = d -> (d.nm))
    g_data = filter(d -> d.nm == nm, results_table)
    println(g_data)

    #Select relevant g (5 --> δh = 0.001)
    g_plus = g_data[index].g_plus 
    g_min = g_data[index].g_min

    #Select relevant c (index = 5 --> δh = 0.001)
    c_plus = g_data[index].c_plus 
    c_min = g_data[index].c_min

    #Construction of H_0
    @suppress begin
        global H_01, H_0h = H_parts(0, 1, 1, no, nm, nf)
    end

    #Construction of H (Z, X, R)
    @suppress begin
        if [Z,X,R] == [0,1,1]
            global H_1, H_h = H_01, H_0h
        else
            global H_1, H_h = H_parts(Z, X, R, no, nm, nf)
        end 
    end

    #Selection of the groundstate
    state_0 = data_nm[1].vec

    #Calcluate groundstate energy
    @suppress begin
        global E_0_p = (state_0' * H_01 * state_0) + (h_c + δh) * (state_0' * H_0h * state_0)
        global E_0_m = (state_0' * H_01 * state_0) + (h_c - δh)* (state_0' * H_0h * state_0)
    end

    for nr in 1:length(l2)
        l = l2[nr]
        cas = c2[nr] 
        i = I[nr]
        #Selection of the state (l2, c2, Z, X, R)
        states = select_states(l,cas,Z, data_nm)
        state = states[i]
        #Finding the scaling dimension
        vec = state.vec
        Δ = state.Δ
        println("Δ=",Δ)
        #Calculating the energy at h ± δh
        @suppress begin
            global E_p = (vec' * H_1 * vec) + (h_c + δh)* (vec' * H_h * vec) 
            global E_m = (vec' * H_1 * vec) + (h_c - δh)* (vec' * H_h * vec)
        end 
        f_ope =  (E_p - E_0_p - (E_m - E_0_m) - (c_plus - c_min)* Δ)/(g_plus-g_min)
        println("f=",f_ope,"\n Z=",Z,"X=",X,"R=",R,"c2=",cas,"l2=",l,"nm=",nm)
        push!(OPEs, (
            nm    = nm,
            Z     = Z,
            X     = X,
            R     = R,
            l2    = l,
            c2    = cas,
            index = i,
            δh    = δh,
            h_c   = h_c,
            Δ     = Δ,
            f_ope = f_ope,
        ))
    end 
end
println(OPEs)

# ── Pretty print ──────────────────────────────────────────────────────────────
println("─"^80)
println("nm    Z    X    R    l2    c2    δh          Δ           f_ope")
println("─"^80)

for i in 1:length(OPEs)
    local r = OPEs[i]
    @printf "%-4d  %-3d  %-3d  %-3d  %-4d  %-4d  %-8.4f  %-10.6f  %-12.6f\n" r.nm r.Z r.X r.R r.l2 r.c2 r.δh r.Δ r.f_ope
    println("─"^80)
end

#Opslaan
outdir  = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/julia/OPE/results"
isdir(outdir) || mkpath(outdir)
outfile = joinpath(outdir, file)

@save outfile OPEs
println("Opgeslagen in: ", outfile)
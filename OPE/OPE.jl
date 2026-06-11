using JLD2
using Dates
using FuzzifiED
using LinearAlgebra
using Printf
using Glob
using SparseArrays  

include("Functions.jl")

#Data from the states
result_dir = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/HPC/O3"
all_data   = load_all_results(result_dir)
sort!(all_data, by = d -> (d.nm, d.E))

#Data from the g
g_file = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/julia/OPE/results/g_results.jld2"
@load g_file results_table
g_data = results_table
#println(g_data)
sort!(g_data, by = d -> (d.nm))

# System size
nm = 7
no = 4*nm
nf = 4
h_c = 14.992 + 6.59 * nm ^ (-2.16 / 2)
R = (nm)^(1/2)

#Operator we want to study
Z=1
X=1
R=1
l2 = 0
c2 = 2
I = 1 #order in the symmetry sector

#Select relevant data, select ground state
data_nm = filter(d -> d.nm == nm, all_data)


#Find vectors of σ and ∂σ
state_σ = select_states(0,2,1,1,1, data_nm) 
state_∂σ = select_states(2,2,1,1,-1, data_nm)

σ  = first(sort(state_σ,  by = d -> d.Δ))
∂σ = first(sort(state_∂σ, by = d -> d.Δ))

eig_σ_sec  = σ.vec
eig_∂σ_sec = ∂σ.vec


#Select g and δh
δh = δh_nm6 = [r.δh for r in g_data if r.nm == nm]
g_plus = [r.g_plus for r in g_data if r.nm == nm]
g_min = [r.g_min for r in g_data if r.nm == nm]
println(g_min)

#Build Hamiltonian for O

H_O1, H_Oh = H_parts(Z, X, R, no, nm, nf)
states_O = select_states(l2,c2,Z,X,R, data_nm) 
state_O = sort(states_O, by = d -> d.Δ)[I]
vec_O = state_O.vec
println(state_O.Δ)

#Build Hamiltonian for vacuum
H_01, H_0h = H_parts(0, 1, 1, no, nm, nf)
eig_0_sec = data_nm[1].vec

results_table = []

for i in 1:length(δh)
    delta_h = δh[i]
    delta_g = g_plus[i] - g_min[i]

    E_plus = (vec_O' * H_O1 * vec_O) + (h_c + delta_h) * (vec_O' * H_Oh * vec_O)
    E_min = (vec_O' * H_O1 * vec_O) + (h_c - delta_h) * (vec_O' * H_Oh * vec_O)

    E_plus_0 = real(eig_0_sec' * H_01 * eig_0_sec) + (h_c + delta_h) * real(eig_0_sec' * H_0h * eig_0_sec)
    E_min_0  = real(eig_0_sec' * H_01 * eig_0_sec) + (h_c - delta_h) * real(eig_0_sec' * H_0h * eig_0_sec)

    f = (E_plus - E_plus_0 - E_min + E_min_0) / delta_g

    push!(results_table, (
        nm     = nm,
        δh     = delta_h,
        h_c    = h_c,
        δ_g = delta_g,
        E_plus = E_plus,
        E_min = E_min,
        f = f,
    ))
end 

println("─"^75)
@printf "%-4s  %-8s  %-10s  %-12s  %-12s %-12s %-12s\n" "nm" "δh" "h_c" "δg" "E_plus" "E_min" "f"
println("─"^75)

current_nm = -1
for r in results_table
    if r.nm != current_nm
        global current_nm = r.nm
        println("─"^75)
    end
    @printf "%-4d  %-8.4f  %-10.4f  %-12.6f  %-12.6f %-12.6f %-12.6f\n" r.nm r.δh r.h_c r.δ_g r.E_plus r.E_min r.f
end
println("─"^75)






#=


c, gr = solve_c_gε_fast(h_c, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0_sec, eig_σ_sec, eig_∂σ_sec)

println(data_nm[1].E)



eig_0_sec = data_nm[1].vec

H_01, H_0h = H_parts(0, 1, 1)
H_σ1, H_σh = H_parts(1, 1, 1)
H_∂σ1, H_∂σh = H_parts(1, 1, -1)


state_ϵ = select_states(0,0,0,1,1, data_nm) 

ϵ = sort(state_ϵ, by = d -> d.Δ)[2]

eig_ϵ_sec = ϵ.vec

E_0_ref = data_nm[1].E

println(E_0_ref)

E_plus = (eig_0_sec' * H_01 * eig_0_sec) + (h_c+δh) * (eig_0_sec' * H_0h * eig_0_sec)
E_min = (eig_0_sec' * H_01 * eig_0_sec) + (h_c-δh) * (eig_0_sec' * H_0h * eig_0_sec)

E_plus_ϵ = (eig_ϵ_sec' * H_01 * eig_ϵ_sec) + (h_c+δh) * (eig_ϵ_sec' * H_0h * eig_ϵ_sec)
E_min_ϵ = (eig_ϵ_sec' * H_01 * eig_ϵ_sec) + (h_c-δh) * (eig_ϵ_sec' * H_0h * eig_ϵ_sec)

δE_plus = E_plus_ϵ - E_0_ref
δE_min = E_min_ϵ - E_0_ref

f_εεε = (δE_plus - δE_min) / (g_plus - g_min) 
println("f_εεε = ", f_εεε)

println("g_plus = ", g_plus) 
println("g_min = ", g_min)
println("g =", g)

println("ϵ: Δ=", ϵ.Δ)  # moet ≈ 1.56 zijn

println("h_c:    ", h_c)
println("g bij h_c: ", g)  # moet ≈ 0 zijn


# FOUT: gebruikt H_01 (Z=0 sector) voor ϵ (ook Z=0, klopt)
# maar enrg_0 moet de ABSOLUTE grondtoestand zijn (Z=0, laagste)
# controleer of data_nm[1] echt het vacuum is:
println("vacuum: Z=", data_nm[1].Z, " X=", data_nm[1].X, " R=", data_nm[1].R, " E=", data_nm[1].E)
println("ϵ:      Z=", ϵ.Z,          " X=", ϵ.X,          " R=", ϵ.R,          " E=", ϵ.E)

# En: δE_plus moet berekend worden bij h_c + δh
# met g_plus ook bij h_c + δh
# dus de afgeleide is correct MAAR alleen als g_ε lineair is in h rond h_c

# Betere aanpak: sweep over meerdere h en fit de helling
h_sweep = range(h_c - 0.3, h_c + 0.3, length=10)
δE_vals = Float64[]
gε_vals = Float64[]

println("c/R = ", c)
println("δE_ϵ = ", real(ϵ.E - data_nm[1].E))
println("(c/R)*Δ_ϵ = ", c * ϵ.Δ)
println("g*f = ", g * f_εεε)
println("check: (c/R)*Δ_ϵ + g*f = ", c * ϵ.Δ + g * f_εεε)
println("moet zijn:               ", real(ϵ.E - data_nm[1].E))

println("f_εεε berekend: ", f_εεε)        # 1.10
println("f_εεε verwacht: ", 0.53)          # bootstrap
println("verhouding: ", f_εεε / 0.53)      # ≈ 2?

# Controleer ook f_σεσ
# Die hebben we als input gebruikt: f_σεσ = 0.525
# Maar wat krijgen we als output voor σ zelf?
E_0_hc  = real(eig_0_sec' * H_01 * eig_0_sec)  + h_c * real(eig_0_sec' * H_0h * eig_0_sec)
E_σ_hc  = real(eig_σ_sec' * H_σ1 * eig_σ_sec)  + h_c * real(eig_σ_sec' * H_σh * eig_σ_sec)
δE_σ_hc = real(E_σ_hc - E_0_hc)

f_σεσ_check = (δE_σ_hc - c * σ.Δ) / g
println("f_σεσ input:    0.525")
println("f_σεσ berekend: ", f_σεσ_check)

# Wat is de laagste energie in sector (Z=1, X=1, R=1)?
data_Z1 = filter(d -> d.Z==1 && d.X==1 && d.R==1, data_nm)
E_0_Z1 = minimum(d.E for d in data_Z1)
println("laagste E in Z=1 sector: ", E_0_Z1)
println("σ.E: ", σ.E)
println("δE_σ correct: ", σ.E - data_nm[1].E)

# Voor ε: vacuum is ook Z=0
println("δE_ϵ / 2 = ", real(ϵ.E - data_nm[1].E) / 2)
println("f_εεε / 2 = ", f_εεε / 2)

println("referentie in solve_c_gε: E_0 = ", 
    real(eig_0_sec' * H_01 * eig_0_sec) + h_c * real(eig_0_sec' * H_0h * eig_0_sec))
println("data_nm[1].E = ", data_nm[1].E)
println("zijn gelijk: ", abs(real(eig_0_sec' * H_01 * eig_0_sec) + 
    h_c * real(eig_0_sec' * H_0h * eig_0_sec) - data_nm[1].E) < 1e-6)

=#
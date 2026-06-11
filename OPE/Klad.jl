#Needed: Functions.jl
#Needed: Data for states and date for g

using JLD2
using Dates
using FuzzifiED
using LinearAlgebra
using Printf
using Glob
using SparseArrays
using Suppressor

#Functions
include("Functions.jl")

#Data
result_dir = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/HPC/O3"
all_data   = load_all_results(result_dir)
sort!(all_data, by = d -> (d.nm, d.E))

nm = 7
data_nm = filter(d -> d.nm == nm, all_data)

l2 =0
c2 =2
Z = 1

state = filter(d ->
    isapprox(d.l2, l2, atol=0.4) &&
    isapprox(d.c2, c2, atol=0.4) && 
    d.Z==Z,
    data_nm)





#=
using JLD2
using Dates
using FuzzifiED
using LinearAlgebra
using Printf
using Glob
using SparseArrays

# System size
nm = 5
no = 4*nm
nf = 4
h_c = 14.992 + 6.59 * nm ^ (-2.16 / 2)
R = (nm)^(1/2)



include("Functions.jl")

#Data from the states
result_dir = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/HPC/O3"
all_data   = load_all_results(result_dir)
sort!(all_data, by = d -> (d.nm, d.E))


data_nm = filter(d -> d.nm == nm, all_data)

#Data from the g
g_file = "C:/Users/marie/OneDrive/Bureaublad/UGent/thesis/julia/OPE/results/g_results.jld2"
@load g_file results_table
g_data = results_table
#println(g_data)
sort!(g_data, by = d -> (d.nm))



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

 
c, g = solve_c_gε_fast(h_c , H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0_sec, eig_σ_sec, eig_∂σ_sec)

c_plus, g_plus = solve_c_gε_fast(h_c + 0.005, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0_sec, eig_σ_sec, eig_∂σ_sec)
c_min, g_min = solve_c_gε_fast(h_c -0.005, H_01, H_0h, H_σ1, H_σh, H_∂σ1, H_∂σh, eig_0_sec, eig_σ_sec, eig_∂σ_sec)

E_0_p = (eig_0_sec' * H_01 * eig_0_sec) + (h_c + 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_σ_p = (eig_σ_sec' * H_σ1 * eig_σ_sec) + (h_c + 0.005)* (eig_σ_sec' * H_σh * eig_σ_sec) 
E_0_m = (eig_0_sec' * H_01 * eig_0_sec) + (h_c - 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_σ_m = (eig_σ_sec' * H_σ1 * eig_σ_sec) + (h_c - 0.005)* (eig_σ_sec' * H_σh * eig_σ_sec) 

println(g, "  ", g_plus, "  ", g_min)

f = (E_σ_p - E_0_p - (E_σ_m - E_0_m))/(g_plus-g_min)


println(f)

f_nieuw =  (E_σ_p - E_0_p - (E_σ_m - E_0_m) - (c_plus - c_min)*0.518936)/(g_plus-g_min)
println("op hoop van zegen", f_nieuw)

println(σ.Δ)

println(∂σ.Δ)

state_ϵ = select_states(0,0,0,1,1, data_nm)
println(state_ϵ[2].Δ) #ϵ

ϵ_vect = state_ϵ[2].vec

E_0_p = (eig_0_sec' * H_01 * eig_0_sec) + (h_c + 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_ϵ_p = (ϵ_vect' * H_01 * ϵ_vect) + (h_c + 0.005)* (ϵ_vect' * H_0h * ϵ_vect) 
E_0_m = (eig_0_sec' * H_01 * eig_0_sec) + (h_c - 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_ϵ_m = (ϵ_vect' * H_01 * ϵ_vect) + (h_c - 0.005)* (ϵ_vect' * H_0h * ϵ_vect) 


f_nieuw =  (E_ϵ_p - E_0_p - (E_ϵ_m - E_0_m) - (c_plus - c_min)*1.59488)/(g_plus-g_min)
println("op hoop van zegen", f_nieuw)
println("g =",g)


########
#Nieuwe primary t(2)
H_t1, H_th = H_parts(0, 1, 1, no, nm, nf)
state_t = select_states(0,6,0,1,1, data_nm)

t_vect = state_t[1].vec

E_0_p = (eig_0_sec' * H_01 * eig_0_sec) + (h_c + 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_t_p = (t_vect' * H_t1 * t_vect) + (h_c + 0.005)* (t_vect' * H_th * t_vect) 
E_0_m = (eig_0_sec' * H_01 * eig_0_sec) + (h_c - 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_t_m = (t_vect' * H_t1 * t_vect) + (h_c - 0.005)* (t_vect' * H_th * t_vect)
f_nieuw =  (E_t_p - E_0_p - (E_t_m - E_0_m) - (c_plus - c_min)*1.20954)/(g_plus-g_min)
println("op hoop van zegen, t(2) = ", f_nieuw)
println("g =",g)


########
#Nieuwe primary j
H_j1, H_jh = H_parts(0, -1, -1, no, nm, nf)
state_j = select_states(2,2,0,-1,-1, data_nm)

j_vect = state_j[1].vec

E_0_p = (eig_0_sec' * H_01 * eig_0_sec) + (h_c + 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_j_p = (j_vect' * H_j1 * j_vect) + (h_c + 0.005)* (j_vect' * H_jh * j_vect) 
E_0_m = (eig_0_sec' * H_01 * eig_0_sec) + (h_c - 0.005)* (eig_0_sec' * H_0h * eig_0_sec)
E_j_m = (j_vect' * H_j1 * j_vect) + (h_c - 0.005)* (j_vect' * H_jh * j_vect)
f_nieuw =  (E_j_p - E_0_p - (E_j_m - E_0_m) - (c_plus - c_min)*2)/(g_plus-g_min)
println("op hoop van zegen, j = ", f_nieuw)
println("g =",g)
=#
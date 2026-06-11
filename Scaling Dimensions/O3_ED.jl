# This example calculates the spectrum of O(3) Wilson-Fisher CFT.
# This example takes the model from 2510.09755 
# and reproduces partly Tables I and II, Figures 1 and 2.
# On my table computer, this calculation takes 6.184 s
using JLD2
using Dates
using FuzzifiED
using LinearAlgebra
FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)
outdir = joinpath(pwd(), "results")
#nm = parse(Int, ARGS[1])

nm = 6 # Number of magnetic flux quanta ~ system size
nf = 4 # Number of 'flavours' (internal states)
no = nm * nf # Total number of orbitals 


# Only one non-zero element in the 4x4 matrix at (4,4)
# Flavour 4

mat_0 = diagm([0, 0, 0, 1.0])

# Adjoint representation 
# Spin-1 representation
# Flavours 1-3
# Only work on 3 dimensions → last row and colum are 0
# Remaining 3x3 part → SO(3) generators ([L_i,L_j] = iϵ_ijk L_k)

mat_A = [
    [ 0  1 0 0 ; 1 0  1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2, # Ax
    [ 0 -1 0 0 ; 1 0 -1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2 * im, # Ay
    [ 1  0 0 0 ; 0 0  0 0 ; 0 0 -1 0 ; 0 0 0 0 ] # Az
]

# Vectro representation 
# O(N) model → interactions between a vectorfield φ
# Couples singlet (flavour 4) to triplet (flavours 1-3)


mat_V = [
    [ 0 0 0  1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0 -1 0 ] / √2, # Vx
    [ 0 0 0 -1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0  1 0 ] / √2 * im, # Vy
    [ 0 0 0  0 ; 0 0 0 1 ; 0 0 0  0 ; 0 1  0 0 ] # Vz
]

# Diagonal Quantum Numbers 
# Particle number
# Angular momentum Lz
# Flavor 1 → charge +1, Flavor 2 and 4 → charge 0, Flavour 3 → charge -1
# Flavors 1-3 get charge +1, Flavor 4 remains 0 (even / odd Parity)

qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, Dict([1 => 1, 3 => -1])), 
    GetFlavQNDiag(nm, nf, Dict([f => 1 for f = 1 : nf - 1]), 0, 2)
]

# Off-diagonal Quantum Numbers
# Rotation
# Flavor parity 

qnf = [
    GetRotyQNOffd(nm, nf),
    GetFlavPermQNOffd(nm, nf, Dict([1 => 3, 3 => 1]))
] ;

# Construction of the Hamiltonian

cfs = Dict{Int64, Confs}()
for Z = 0 : 1
    cfs[Z] = Confs(no, [nm, 0, 0, Z], qnd)
end 

# Value for the critical point

h = 14.992 + 6.59 * nm ^ (-2.16 / 2)

# Physical interaction between particles

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, [6.5, 1.0])
    - 1.4 * 2 * GetDenIntTerms(nm, nf, [6.5, 1.0], mat_V)
    - 2 * h * GetPolTerms(nm, nf, mat_0)
)

# Labels for the states (L^2, Casimir Operator)

tms_l2 = GetL2Terms(nm, nf)
tms_c2 = 4 * GetC2Terms(nm, nf, mat_A) ;
 

results = Dict()
eigenvectors = Dict()
basis_info = Dict()
result = []

for Z in [0, 1], X in [1, -1], R in [1, -1]
    key = (Z, X, R)

    bs = Basis(cfs[Z], [R, X], qnf)
    
    hmt = Operator(bs, tms_hmt)
    hmt_mat = OpMat(hmt)
    enrg, st = GetEigensystem(hmt_mat, 20) #Number of states saved per sector

    basis_info[key] = (
        Z = Z,
        X = X,
        R = R,
        ne = nm,
    )

    l2 = Operator(bs, tms_l2)
    l2_mat = OpMat(l2)
    l2_val = [ st[:, i]' * l2_mat * st[:, i] for i in eachindex(enrg)]

    c2 = Operator(bs, tms_c2)
    c2_mat = OpMat(c2)
    c2_val = [ st[:, i]' * c2_mat * st[:, i] for i in eachindex(enrg)]

    results[key] = (
        energy = enrg,
        l2 = l2_val,
        c2 = c2_val
    )

    eigenvectors[key] = st

    for i in eachindex(enrg)
        push!(result, [enrg[i], l2_val[i], c2_val[i], X, R, Z])
    end

    jobid     = get(ENV, "PBS_JOBID", "local")
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    tmpfile   = joinpath(outdir, "O3_nm$(nm)_$(jobid)_checkpoint.jld2")
    @save tmpfile nm nf no h result results basis_info
    println("Checkpoint saved after (Z=$Z, X=$X, R=$R) at $timestamp")
end


sort!(result, by = st -> real(st[1]))
enrg_0 = result[1][1]
enrg_ϕ = filter(st -> st[2] ≈ 0 && st[3] ≈ 2 && st[4] ≈ 1, result)[1][1]
spec = [ round.([ 0.518936 * (st[1] - enrg_0) / (enrg_ϕ - enrg_0) ; st] .+ √eps(Float64), digits = 6) for st in result ]
display(permutedims(hcat(spec...)))


# Scaling dimension → Energy → Angular momentum → Casimir → Parity

outdir = joinpath(pwd(), "results")
isdir(outdir) || mkpath(outdir)

jobid = get(ENV, "PBS_JOBID", "local")
timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
outfile = joinpath(outdir, "O3_nm$(nm)_$(jobid)_$(timestamp).jld2")
@save outfile nm nf no h enrg_0 enrg_ϕ spec result results eigenvectors basis_info
#@save outfile nm nf no h enrg_0 enrg_ϕ spec result results basis_info
println("Data saved to: ", outfile)
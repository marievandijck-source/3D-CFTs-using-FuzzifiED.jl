using FuzzifiED
using Dates
using LinearAlgebra
using ITensors, ITensorMPS
using Printf
using JLD2


t_start = now()

# BLAS parallel (matrix operaties), Julia sequentieel (DMRG loops)
BLAS.set_num_threads(8)
ITensors.Strided.set_num_threads(1)
ITensors.BLAS.set_num_threads(8)
println("BLAS threads: ", BLAS.get_num_threads())
println("Julia threads: ", Threads.nthreads()) 


FuzzifiED.ElementType = Float64
≈(x, y) = abs(x - y) < √eps(Float64)

#=
nm = parse(Int, ARGS[1])
even = parse(Int, ARGS[2])
odd = parse(Int, ARGS[3])
=#

nm = 14
even = 2
odd = 1

nf = 4
no = nm * nf


function get_dmrg_kwargs(nm)
    if nm <= 10
        return (nsweeps=10,
                maxdim=[10, 20, 50, 100, 200, 200, 300, 300, 300, 300],
                noise= [1e-3,3e-4,1e-5,3e-6,1e-6,  0,   0,   0,   0,   0],
                cutoff=[1e-9],
                eigsolve_tol=1e-5,
                eigsolve_maxiter=8)
    elseif nm <= 14
        return (nsweeps=12,
                maxdim=[10,20,50,100,200,400,600,600,600,600,600,600],
                noise=[1e-4,3e-5,1e-5,3e-6,1e-6,1e-7,1e-8,1e-9,0,0,0,0],
                cutoff=[1e-10],
                eigsolve_tol=1e-6,
                eigsolve_maxiter=10)
    elseif nm <= 15
        return (nsweeps=14,
                maxdim=[10,20,50,100,200,400,800,900,900,900,900,900,900,900],
                noise=[1e-4,3e-5,1e-5,3e-6,1e-6,1e-7,1e-8,1e-9,0,0,0,0,0,0],
                cutoff=[1e-11],
                eigsolve_tol=1e-7,
                eigsolve_maxiter=10)
    else
        return (nsweeps=16,
                maxdim=[10,20,50,100,200,400,800,1200,1500,2000,2000,2000,2000,2000,2000,2000],
                noise=[1e-4,3e-5,1e-5,3e-6,1e-6,1e-7,1e-8,1e-9,1e-10,0,0,0,0,0,0,0],
                cutoff=[1e-12],
                eigsolve_tol=1e-8,
                eigsolve_maxiter=10)
    end
end

# ── DMRG instellingen ─────────────────────────────────────────
dmrg_kwargs = get_dmrg_kwargs(nm)
weight = 100.0

# ── Representatiematrices ─────────────────────────────────────
mat_V = [
    [ 0 0 0  1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0 -1 0 ] / √2,
    [ 0 0 0 -1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0  1 0 ] / √2 * im,
    [ 0 0 0  0 ; 0 0 0 1 ; 0 0 0  0 ; 0 1  0 0 ]
]
mat_A = [
    [ 0  1 0 0 ; 1 0  1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2,
    [ 0 -1 0 0 ; 1 0 -1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2 * im,
    [ 1  0 0 0 ; 0 0  0 0 ; 0 0 -1 0 ; 0 0 0 0 ]
]
mat_0 = diagm([0, 0, 0, 1.0])

# ── Sites en kwantumgetallen ──────────────────────────────────
sites = GetSites([
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, Dict([1 => 1, 3 => -1])),
    GetFlavQNDiag(nm, nf, Dict([f => 1 for f = 1 : nf - 1]), 0, 2)
])

# ── Hamiltoniaan ──────────────────────────────────────────────
h = 14.992 + 6.59 * nm ^ (-2.16 / 2)

tms_hmt = SimplifyTerms(
    GetDenIntTerms(nm, nf, [6.5, 1.0])
    - 1.4 * 2 * GetDenIntTerms(nm, nf, [6.5, 1.0], mat_V)
    - 2 * h * GetPolTerms(nm, nf, mat_0)
)
hmt = MPO(OpSum(tms_hmt), sites)

# ── Observables ───────────────────────────────────────────────
l2 = MPO(OpSum(GetL2Terms(nm, nf)),          sites)
c2 = MPO(OpSum(4 * GetC2Terms(nm, nf, mat_A)), sites)

# ── Hulpfuncties voor begintoestanden ─────────────────────────
f_of(o)    = (o - 1) % nf + 1
m_of(o)    = (o - 1) ÷ nf + 1
base_f(m)  = isodd(m) ? 1 : 3
alt_f(m)   = isodd(m) ? 3 : 1
site(m, f) = (m - 1) * nf + f

# Basisbezetting: Ne=nm, Lz=0, Q=0, Z=0
# f=1 op m=1..nm÷2, f=3 op m=nm÷2+1..nm
cfi_0 = [ (f_of(o)==1 && m_of(o) <= nm÷2) ? 1 :
          (f_of(o)==3 && m_of(o) >  nm÷2) ? 1 : 0
          for o = 1:no ]

cfi_phi = copy(cfi_0)
cfi_phi[site(nm÷2 + 1, 3)] = 0   # verwijder één f=3
cfi_phi[site(nm÷2 + 1, 4)] = 1   # vervang door f=4


sti_test = MPS(sites, string.(cfi_0))
println("Z van cfi_0: ", val(flux(sti_test)[4]))

sti_test_phi = MPS(sites, string.(cfi_phi))
println("Z van cfi_phi: ", val(flux(sti_test_phi)[4]))
#=
# φz sector: wissel f=base(m=1) → alt_f(m=1)
cfi_z = copy(cfi_0)
cfi_z[site(1, base_f(1))] = 0
cfi_z[site(1, alt_f(1))]  = 1

# φ± sector: wissel f=base(m=1) → f=2
cfi_pm = copy(cfi_0)
cfi_pm[site(1, base_f(1))] = 0
cfi_pm[site(1, 2)]          = 1

# J sector: Lz≠0
cfi_J = copy(cfi_0)
cfi_J[site(1, 1)] = 0
cfi_J[site(1, 3)] = 1
cfi_J[site(3, 1)] = 0
cfi_J[site(4, 1)] = 1

# Tμν sector
cfi_Tmunu = copy(cfi_0)
cfi_Tmunu[site(3, 1)] = 0
cfi_Tmunu[site(5, 1)] = 1

# ∂²φ sector
cfi_delphi2 = copy(cfi_z)
cfi_delphi2[site(3, 1)] = 0
cfi_delphi2[site(5, 1)] = 1
=#




# ── Begintoestanden als MPS ───────────────────────────────────
sti_0 = MPS(sites, string.(cfi_0))

#sti_z = MPS(sites, string.(cfi_z))
sti_phi = MPS(sites, string.(cfi_phi))
#sti_pm = MPS(sites, string.(cfi_pm))
#sti_J = MPS(sites, string.(cfi_J))
#sti_Tmunu = MPS(sites, string.(cfi_Tmunu))
#sti_delphi2 = MPS(sites, string.(cfi_delphi2))

# ── Begintoestanden en aantal toestanden per flux-sector ──────
# Elke begintoestand heeft een andere flux (QN-combinatie)
# n_states: hoeveel aangeslagen toestanden per sector
init_list = [
    (sti_0,       even),
    (sti_phi,       odd),
]

println("flux sti_0: ", flux(sti_0))
println("flux sti_phi", flux(sti_phi))
#println("flux sti_J:  ", flux(sti_J))
#println("flux sti_Tmunu: ", flux(sti_Tmunu))
#println("flux sti_delphi2: ", flux(sti_delphi2))


# ── DMRG lus ─────────────────────────────────────────────────
all_states = []

all_mps = []

for (sti, n_states) in init_list
    println("\n===== flux = $(flux(sti)) =====")

    prev = Vector{MPS}()

    for i in 1:n_states
        println("  Toestand $i/$n_states")

        E, st = isempty(prev) ?
            dmrg(hmt, sti; dmrg_kwargs...) :
            dmrg(hmt, prev, sti; weight=weight, dmrg_kwargs...)

        push!(prev, st)

        l2_val = real(inner(st', l2, st))
        c2_val = real(inner(st', c2, st))

        push!(all_states, (
            E      = E,
            l2     = l2_val,
            c2     = c2_val,
            flux   = flux(st),
            Z = Int(val(flux(st)[4])),
        ))

        push!(all_mps, st)

        @printf "    E=%.6f  L²=%.2f  C₂=%.2f\n" E l2_val c2_val
    end
end

perm = sortperm(all_states, by = st -> st.E)
all_states = all_states[perm]
all_mps    = all_mps[perm]

gs       = all_states[1]
gs_mps   = all_mps[1]
gs_file  = "O3_dmrg_nm$(nm)_gs.jld2"
@save gs_file nm gs gs_mps
println("Grondtoestand opgeslagen in $gs_file")

# ── Sorteer op energie ────────────────────────────────────────
sort!(all_states, by = st -> st.E)

# ── Referentie-energieën ──────────────────────────────────────
enrg_0 = all_states[1].E

# φ toestand: C₂ ≈ 2, L² ≈ 0
phi_states = filter(st -> abs(st.c2 - 2) < 0.5 && abs(st.l2) < 0.5, all_states)
if isempty(phi_states)
    @warn "Geen φ-toestand gevonden met C₂≈2, L²≈0 — check je sectoren"
    enrg_ϕ = all_states[2].E  # fallback
else
    enrg_ϕ = phi_states[1].E
end

# ── Print spectrum ────────────────────────────────────────────
println("\n===== Spectrum =====")
@printf "%-8s  %-12s  %-8s  %-8s %-8s\n" "Δ" "E" "L²" "C₂" "Z"
println("─"^45)
for st in all_states
    Δ = 0.518936 * (st.E - enrg_0) / (enrg_ϕ - enrg_0)
    @printf "%-8.4f  %-12.6f  %-8.3f  %-8.3f %-d\n" Δ st.E st.l2 st.c2 st.Z
end

# ── Opslaan ───────────────────────────────────────────────────
timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
outfile   = "O3_dmrg_nm$(nm)_$(timestamp).jld2"
@save outfile nm nf no h enrg_0 enrg_ϕ all_states
println("\nOpgeslagen in $outfile")



t_end = now()
elapsed = t_end - t_start
println("Totale rekentijd: ", canonicalize(elapsed))
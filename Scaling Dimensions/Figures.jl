using JLD2
using Plots
using Plots.PlotMeasures
using Pkg
using LaTeXStrings
using Colors
using ColorSchemes


dmrg_files = filter(f -> startswith(f, "O3_dmrg_nm") &&
                         endswith(f, ".jld2"),
                    readdir(@__DIR__))

dmrg_data = []

for f in dmrg_files

    local nm, all_states, enrg_0, enrg_ϕ

    @load joinpath(@__DIR__, f) nm all_states enrg_0 enrg_ϕ

    for st in all_states

        Δ = 0.518936 * (st.E - enrg_0) / (enrg_ϕ - enrg_0)

        push!(dmrg_data, (
            nm = nm,
            Δ  = Δ,
            E  = st.E,
            l2 = st.l2,
            c2 = st.c2,
            Z  = st.Z,
            method = :DMRG,
        ))
    end
end

println("Geladen DMRG toestanden: ", length(dmrg_data))

# ── Parameters ────────────────────────────────────────────────
ω = 0.7668
x_of_nm(nm) = sqrt(Float64(nm))^(-ω)

# ── Data laden ────────────────────────────────────────────────
@load joinpath(@__DIR__, "O3_combined2.jld2") all_data
println("Geladen: ", length(all_data), " toestanden")

nms_ed   = unique([d.nm for d in all_data])
nms_dmrg = unique([d.nm for d in dmrg_data])
nms      = sort(unique(vcat(nms_ed, nms_dmrg)))

all_x = x_of_nm.(nms)
xmin  = 0
xmax  = 0.6

# ── Bootstrap waarden ─────────────────────────────────────────
# (L, c2, Z_even, Δ, label)
bootstrap = [
    (L=0, c2=2, Z_even=false, Δ=0.518936, label=L"\sigma"),
    (L=0, c2=0, Z_even=true,  Δ=1.59488,  label=L"\epsilon"),
    (L=0, c2=0, Z_even=true,  Δ=3.76680,  label=L"\epsilon'"),
    (L=0, c2=6, Z_even=true,  Δ=1.20954,  label=L"t_{(2)}"),
  

]

bootstrap_dc = [
    (L=1, c2=2, Z_even=false, Δ=1.518936, label=L"\partial\sigma"),
    (L=0, c2=2, Z_even=false, Δ=2.518936, label=L"\partial\partial\sigma^{(0)}"),
    (L=2, c2=2, Z_even=false, Δ=2.518936, label=L"\partial\partial\sigma^{(2)}"),
    (L=1, c2=2, Z_even=false, Δ=3.518936, label=L"\partial\partial\partial\sigma^{(1)}"),
    (L=0, c2=2, Z_even=false, Δ=4.518936, label=L"\partial\partial\partial\partial\sigma^{(0)}"),
    (L=2, c2=2, Z_even=false, Δ=4.518936, label=L"\partial\partial\partial\partial\sigma^{(2)}"),

    (L=1, c2=0, Z_even=true,  Δ=2.59488,  label=L"\partial\epsilon"),
    (L=2, c2=0, Z_even=true,  Δ=3.59488,  label=L"\partial\partial\epsilon^{(2)}"),
    (L=0, c2=0, Z_even=true,  Δ=3.59488,  label=L"\partial\partial\epsilon^{(0)}"),
    (L=1, c2=0, Z_even=true,  Δ=4.59488,  label=L"\partial\partial\partial\epsilon^{(1)}"),
    (L=1, c2=0, Z_even=true,  Δ=4.76680,  label=L"\partial\epsilon'"),

    (L=0, c2=6, Z_even=true,  Δ=3.20954,  label=L"\partial\partial t_{(2)}^{(0)}"),
    (L=2, c2=6, Z_even=true,  Δ=3.20954,  label=L"\partial\partial t_{(2)}^{(2)}"),
    (L=1, c2=6, Z_even=true,  Δ=2.20954,  label=L"\partial\partial t_{(2)}^{(1)}"),
    (L=1, c2=6, Z_even=true,  Δ=4.20954,  label=L"\partial\partial\partial t_{(2)}^{(1)}"),
]

primaries = [
    (L=0, c2=0, Z_even=true,  Δ=0.0,      label=L"\mathbf{1}"),
    (L=2, c2=0, Z_even=true,  Δ=3.0,      label=L"T_{\mu\nu}"),
    (L=1, c2=2, Z_even=true,  Δ=2.0,      label=L"j_\mu"),

]

descendants = [
(L=1, c2=2, Z_even=true,  Δ=3.0,      label=L"\epsilon\partial j_\mu"),
(L=1, c2=2, Z_even=true,  Δ=4.0,      label=""),
(L=2, c2=2, Z_even=true,  Δ=3.0,      label=L"\partial j_\mu"),
(L=2, c2=2, Z_even=true,  Δ=4.0,      label=L"\epsilon\partial j_\mu"),

(L=2, c2=0, Z_even=true,  Δ=4.0,      label=L"\epsilon\partial T_{\mu\nu}"),
]

labels = [
(L=1, c2=0, Z_even=false,  Δ=2.8,      label=L"\epsilon_\mu"),
(L=1, c2=2, Z_even=true,  Δ=3.7,      label=L"\phi_\mu"),
(L=0, c2=2, Z_even=false,  Δ=3.7,      label=L"\sigma'"),
(L=1, c2=2, Z_even=false,  Δ=2.8,      label=L"\sigma_\mu"),
(L=2, c2=2, Z_even=false,  Δ=3.6,      label=L"\sigma_{\mu\nu}"),

(L=1, c2=6, Z_even=false,  Δ=2.68,      label=L"t_{(2)\mu}"),
(L=2, c2=6, Z_even=false,  Δ=3.29,      label=L"t'_{(2)\mu\nu}"),
(L=0, c2=6, Z_even=true,  Δ=3.6,      label=L"t'_{(2)}"),
(L=2, c2=6, Z_even=true,  Δ=2.9,      label=L"t''_{(2)\mu\nu}"),
(L=0, c2=6, Z_even=true,  Δ=4.3,      label=L"t''_{(2)}"),
]

labels_dc = [
(L=1, c2=0, Z_even=false,  Δ=3.8,      label=L"\epsilon\partial\epsilon_\mu"),
(L=0, c2=0, Z_even=false,  Δ=3.8,      label=L"\partial\epsilon_\mu^{(0)}"),
(L=2, c2=0, Z_even=false,  Δ=3.5,      label=L"\partial\epsilon_\mu^{(2)}"),
(L=0, c2=6, Z_even=false,  Δ=3.68,      label=L"\partial t_{(2)\mu}^{(0)}"),
(L=1, c2=6, Z_even=false,  Δ=3.68,      label=L"\epsilon\partial t_{(2)\mu}^{(1)}"),
]

palette = collect(cgrad(
[:maroon, :firebrick1, :yellow, :navajowhite],
length(nms),
categorical=true))

nm_to_color = Dict(nm => palette[i] for (i,nm) in enumerate(nms))

# ── Sectoren om te plotten ────────────────────────────────────
# (c2_target, Z_even, titel, bestandsnaam)
sectors = [
    (0, true,
     L"S = 0^{+}\;(\mathrm{singlet},\ Z_{\mathrm{even}})",
     "O3_S0plus"),

    (0, false,
     L"S = 0^{-}\;(\mathrm{singlet},\ Z_{\mathrm{odd}})",
     "O3_S0min"),

    (2, true,
     L"S = 1^{+}\;(\mathrm{triplet},\ Z_{\mathrm{even}})",
     "O3_S1plus"),

    (2, false,
     L"S = 1^{-}\;(\mathrm{triplet},\ Z_{\mathrm{odd}})",
     "O3_S1min"),

    (6, true,
     L"S = 2^{+}\;(\mathrm{quintet},\ Z_{\mathrm{even}})",
     "O3_S2plus"),

    (6, false,
     L"S = 2^{-}\;(\mathrm{quintet},\ Z_{\mathrm{odd}})",
     "O3_S2min"),
]

L_vals  = [0, 1, 2]
l2_vals = [0, 2, 6]

for (c2_target, Z_even, sector_title, fname) in sectors

    # Maak drie subplots naast elkaar: L=0, L=1, L=2
    plts = []

    for (col, (L, l2_target)) in enumerate(zip(L_vals, l2_vals))

        p = plot(
            xlabel     = L"\mathbf{R^{-\omega}}",
            ylabel     = col == 1 ?  L"\mathbf{\Delta}" : "",
            title = L"L = %$L",
            titlefont  = font(11),
            xlims      = (0, 0.6),
            ylims      = (-0.2, 5.0),
            xticks     = ((0.0,0.2,0.4,0.6), (0.0,0.2,0.4,0.6)),
            yticks     = 0:1:5,
            legend = false,
            grid       = true,
            gridalpha  = 0.2,
            framestyle = :box,
            size       = (350, 400),
        )

        # Bootstrap gestippelde lijnen
        for bs in bootstrap
            bs.L == L && bs.Z_even == Z_even && abs(bs.c2 - c2_target) < 0.5 || continue
            hline!(p, [bs.Δ];
                linestyle = :dash,
                linecolor = :steelblue,
                linewidth = 1.5,
                label     = "",
            )
            annotate!(p, xmax + 0.01, bs.Δ,
                text(bs.label, :steelblue, :left, 11))
        end

        for bs in bootstrap_dc
            bs.L == L && bs.Z_even == Z_even && abs(bs.c2 - c2_target) < 0.5 || continue
            hline!(p, [bs.Δ];
                linestyle = :dash,
                linecolor = :tomato,
                linewidth = 1.5,
                label     = "",
            )
            annotate!(p, xmax + 0.01, bs.Δ,
                text(bs.label, :tomato, :left, 11))
        end

        for bs in descendants
            bs.L == L && bs.Z_even == Z_even && abs(bs.c2 - c2_target) < 0.5 || continue
            hline!(p, [bs.Δ];
                linecolor = :tomato,
                linewidth = 1.5,
                label     = "",
            )
            annotate!(p, xmax + 0.01, bs.Δ,
                text(bs.label, :tomato, :left, 11))
        end

        for bs in primaries
            bs.L == L && bs.Z_even == Z_even && abs(bs.c2 - c2_target) < 0.5 || continue
            hline!(p, [bs.Δ];
                linecolor = :steelblue,
                linewidth = 1.5,
                label     = "",
            )
            annotate!(p, xmax + 0.01, bs.Δ,
                text(bs.label, :steelblue, :left, 11))
        end

        for bs in labels
            bs.L == L && bs.Z_even == Z_even && abs(bs.c2 - c2_target) < 0.5 || continue
            annotate!(p, xmax + 0.01, bs.Δ,
                text(bs.label, :steelblue, :left, 11))
        end

        for bs in labels_dc
            bs.L == L && bs.Z_even == Z_even && abs(bs.c2 - c2_target) < 0.5 || continue
            annotate!(p, xmax + 0.01, bs.Δ,
                text(bs.label, :tomato, :left, 11))
        end

        # Data punten gefilterd op L, C₂ en Z
        sector_data = filter(all_data) do d
            abs(d.l2 - l2_target) < 0.2 + 0.3*L &&
            iseven(d.Z) == Z_even &&
            abs(d.c2 - c2_target) < 0.5
        end
        println("nm waarden in sector_data: ", sort(unique([d.nm for d in sector_data])))

        for d in sector_data
            scatter!(p, [x_of_nm(d.nm)], [d.Δ];
                markercolor = nm_to_color[d.nm],
                markerstrokewidth = 0.3,
                markersize  = 5,
                label       = "",
            )
        end

        dmrg_sector_data = filter(dmrg_data) do d
        abs(d.l2 - l2_target) < 0.2 + 0.3*L &&
        iseven(d.Z) == Z_even &&
        abs(d.c2 - c2_target) < 0.5
        end

        println("DMRG nm waarden in sector_data: ",
        sort(unique([d.nm for d in dmrg_sector_data])))

        for d in dmrg_sector_data
        scatter!(p, [x_of_nm(d.nm)], [d.Δ];
            markercolor = nm_to_color[d.nm],
            markerstrokecolor = :black,
            markerstrokewidth = 0.1,
            markershape = :diamond,
            markersize  = 5,
            label       = "",
        )
        end

        push!(plts, p)
    end
    
# ── Apart legendablok ─────────────────────────────────────────
legend_plot = plot(
    framestyle = :box,
    axis       = false,
    grid       = false,
    legend     = :topleft,
    legendfontsize = 9,
    legendrowheight = 20,
    margin = 0mm,
)

for nm in nms
    scatter!(legend_plot, [NaN], [NaN];
        markercolor = nm_to_color[nm],
        markershape = :rect,
        markersize = 10,
        label = "N = $nm",
    )
    # Lege regel als spacer
    scatter!(legend_plot, [NaN], [NaN];
        markercolor = :white,
        markerstrokewidth = 0,
        markersize  = 0,
        label       = " ",   # spatie als lege label
    )

end
scatter!(legend_plot, [NaN], [NaN];
    markercolor = :gray40,
    markerstrokecolor = :black,
    markershape = :circle,
    markersize = 6,
    label = "ED",
)
scatter!(legend_plot, [NaN], [NaN];
    markercolor = :gray40,
    markerstrokecolor = :black,
    markershape = :diamond,
    markersize = 6,
    label = "DMRG",
)

plot!(legend_plot, [NaN, NaN], [NaN, NaN];
    linecolor = :tomato, linestyle = :solid, linewidth = 1.5,
    label = "descendant (exact)",
)
plot!(legend_plot, [NaN, NaN], [NaN, NaN];
    linecolor = :tomato, linestyle = :dash, linewidth = 0.6,
    label = "descendant (CB)",
)
plot!(legend_plot, [NaN, NaN], [NaN, NaN];
    linecolor = :steelblue, linestyle = :solid, linewidth = 1.5,
    label = "primary (exact)",
)
plot!(legend_plot, [NaN, NaN], [NaN, NaN];
    linecolor = :steelblue, linestyle = :dash, linewidth = 0.6,
    label = "primary (CB)",
)


# ── Combineer drie subplots + legendablok ─────────────────────
fig = plot(
    plts[1], plts[2], plts[3], legend_plot;

    layout = @layout([a{0.27w} b{0.27w} c{0.27w} d{0.22w}]),
    size = (1300, 460),

    plot_title = sector_title,
    plot_titlefontsize = 18,
    plot_titlevspan = 0.09,

    margin = 3mm,
    left_margin = 8mm,
    right_margin = 4mm,
    top_margin = 0mm,
    bottom_margin = 5mm,
)
 

    # Voeg legenda entries toe aan laatste subplot


    # Opslaan
    outfile_pdf = joinpath(@__DIR__, fname * "2.pdf")
    outfile_png = joinpath(@__DIR__, fname * "2.png")
    savefig(fig, outfile_pdf)
    savefig(fig, outfile_png)
    println("Opgeslagen: $fname")
end

println("\nKlaar — 6 figuren aangemaakt.")

using JLD2, Plots, Plots.PlotMeasures
using LaTeXStrings, Colors, ColorSchemes

# ─────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────

ω        = 0.7668
x_of_nm(nm) = Float64(nm)^(-ω / 2)

data_dir = joinpath(@__DIR__, "results/L0S0+")
all_data = []
for f in filter(f -> endswith(f, ".jld2"), readdir(data_dir; join = true))
    append!(all_data, load(f)["OPEs"])
end
println("Loaded: ", length(all_data), " states")

nms = sort(unique([d.nm for d in all_data]))

# ─────────────────────────────────────────────
# Color palette (one color per nm)
# ─────────────────────────────────────────────

palette    = collect(cgrad([:maroon, :firebrick1, :yellow, :navajowhite],
                           length(nms); categorical = true))
nm_to_color = Dict(nm => palette[i] for (i, nm) in enumerate(nms))

xmin, xmax = 0.0, 0.6

# ─────────────────────────────────────────────
# Reference lines from conformal bootstrap
# ─────────────────────────────────────────────

bootstrap = [
    (L=0, c2=2, Z_even=false, Δ=0.518936,  label=L"\sigma"),
    (L=0, c2=0, Z_even=true,  Δ=1.59488,   label=L"\epsilon"),
    (L=0, c2=0, Z_even=true,  Δ=3.76680,   label=L"\epsilon'"),
    (L=0, c2=6, Z_even=true,  Δ=1.20954,   label=L"t_{(2)}"),
]
bootstrap_dc = [
    (L=1, c2=2,  Z_even=false, Δ=1.518936,  label=L"\partial\sigma"),
    (L=0, c2=2,  Z_even=false, Δ=2.518936,  label=L"\partial^2\sigma^{(0)}"),
    (L=2, c2=2,  Z_even=false, Δ=2.518936,  label=L"\partial^2\sigma^{(2)}"),
    (L=1, c2=2,  Z_even=false, Δ=3.518936,  label=L"\partial^3\sigma^{(1)}"),
    (L=0, c2=2,  Z_even=false, Δ=4.518936,  label=L"\partial^4\sigma^{(0)}"),
    (L=2, c2=2,  Z_even=false, Δ=4.518936,  label=L"\partial^4\sigma^{(2)}"),
    (L=1, c2=0,  Z_even=true,  Δ=2.59488,   label=L"\partial\epsilon"),
    (L=2, c2=0,  Z_even=true,  Δ=3.59488,   label=L"\partial^2\epsilon^{(2)}"),
    (L=0, c2=0,  Z_even=true,  Δ=3.59488,   label=L"\partial^2\epsilon^{(0)}"),
    (L=1, c2=0,  Z_even=true,  Δ=4.59488,   label=L"\partial^3\epsilon^{(1)}"),
    (L=1, c2=0,  Z_even=true,  Δ=4.76680,   label=L"\partial\epsilon'"),
    (L=0, c2=6,  Z_even=true,  Δ=3.20954,   label=L"\partial^2 t_{(2)}^{(0)}"),
    (L=2, c2=6,  Z_even=true,  Δ=3.20954,   label=L"\partial^2 t_{(2)}^{(2)}"),
    (L=1, c2=6,  Z_even=true,  Δ=2.20954,   label=L"\partial t_{(2)}^{(1)}"),
    (L=1, c2=6,  Z_even=true,  Δ=4.20954,   label=L"\partial^3 t_{(2)}^{(1)}"),
]
primaries = [
    (L=0, c2=0, Z_even=true,  Δ=0.0,  label=L"\mathbf{1}"),
    (L=2, c2=0, Z_even=true,  Δ=3.0,  label=L"T_{\mu\nu}"),
    (L=1, c2=2, Z_even=true,  Δ=2.0,  label=L"j_\mu"),
]
descendants = [
    (L=1, c2=2, Z_even=true,  Δ=3.0,  label=L"\epsilon\partial j_\mu"),
    (L=1, c2=2, Z_even=true,  Δ=4.0,  label=""),
    (L=2, c2=2, Z_even=true,  Δ=3.0,  label=L"\partial j_\mu"),
    (L=2, c2=2, Z_even=true,  Δ=4.0,  label=L"\epsilon\partial j_\mu"),
    (L=2, c2=0, Z_even=true,  Δ=4.0,  label=L"\epsilon\partial T_{\mu\nu}"),
]
labels = [
    (L=1, c2=0, Z_even=false, Δ=2.8,   label=L"\epsilon_\mu"),
    (L=1, c2=2, Z_even=true,  Δ=3.7,   label=L"\phi_\mu"),
    (L=0, c2=2, Z_even=false, Δ=3.7,   label=L"\sigma'"),
    (L=1, c2=2, Z_even=false, Δ=2.8,   label=L"\sigma_\mu"),
    (L=1, c2=6, Z_even=false, Δ=2.68,  label=L"t_{(2)\mu}"),
    (L=2, c2=6, Z_even=false, Δ=3.29,  label=L"t'_{(2)\mu\nu}"),
    (L=0, c2=6, Z_even=true,  Δ=3.6,   label=L"t'_{(2)}"),
    (L=2, c2=6, Z_even=true,  Δ=2.9,   label=L"t''_{(2)\mu\nu}"),
    (L=0, c2=6, Z_even=true,  Δ=4.3,   label=L"t''_{(2)}"),
]
labels_dc = [
    (L=1, c2=0, Z_even=false, Δ=3.8,   label=L"\epsilon\partial\epsilon_\mu"),
    (L=0, c2=0, Z_even=false, Δ=3.8,   label=L"\partial\epsilon_\mu^{(0)}"),
    (L=2, c2=0, Z_even=false, Δ=3.5,   label=L"\partial\epsilon_\mu^{(2)}"),
    (L=0, c2=6, Z_even=false, Δ=3.68,  label=L"\partial t_{(2)\mu}^{(0)}"),
    (L=1, c2=6, Z_even=false, Δ=3.68,  label=L"\epsilon\partial t_{(2)\mu}^{(1)}"),
]

# Helper: draw reference lines and annotations for one subplot
function add_reference_lines!(p, L, c2_target, Z_even)
    match(b) = b.L == L && b.Z_even == Z_even && abs(b.c2 - c2_target) < 0.5

    for b in filter(match, bootstrap)
        hline!(p, [b.Δ]; linestyle=:dash, linecolor=:steelblue, linewidth=1.5, label="")
        annotate!(p, xmax + 0.01, b.Δ, text(b.label, :steelblue, :left, 11))
    end
    for b in filter(match, bootstrap_dc)
        hline!(p, [b.Δ]; linestyle=:dash, linecolor=:tomato, linewidth=1.5, label="")
        annotate!(p, xmax + 0.01, b.Δ, text(b.label, :tomato, :left, 11))
    end
    for b in filter(match, primaries)
        hline!(p, [b.Δ]; linecolor=:steelblue, linewidth=1.5, label="")
        annotate!(p, xmax + 0.01, b.Δ, text(b.label, :steelblue, :left, 11))
    end
    for b in filter(match, descendants)
        hline!(p, [b.Δ]; linecolor=:tomato, linewidth=1.5, label="")
        annotate!(p, xmax + 0.01, b.Δ, text(b.label, :tomato, :left, 11))
    end
    for b in filter(match, labels)
        annotate!(p, xmax + 0.01, b.Δ, text(b.label, :steelblue, :left, 11))
    end
    for b in filter(match, labels_dc)
        annotate!(p, xmax + 0.01, b.Δ, text(b.label, :tomato, :left, 11))
    end
end

# ─────────────────────────────────────────────
# Sectors to plot: (c2, Z_even, title, filename)
# ─────────────────────────────────────────────

sectors = [
    (0, true,  L"S = 0^{+}\;(\mathrm{singlet},\ Z_\mathrm{even})", "O3_S0plus"),
    (0, false, L"S = 0^{-}\;(\mathrm{singlet},\ Z_\mathrm{odd})",  "O3_S0min"),
    (2, true,  L"S = 1^{+}\;(\mathrm{triplet},\ Z_\mathrm{even})", "O3_S1plus"),
    (2, false, L"S = 1^{-}\;(\mathrm{triplet},\ Z_\mathrm{odd})",  "O3_S1min"),
    (6, true,  L"S = 2^{+}\;(\mathrm{quintet},\ Z_\mathrm{even})", "O3_S2plus"),
    (6, false, L"S = 2^{-}\;(\mathrm{quintet},\ Z_\mathrm{odd})",  "O3_S2min"),
]

L_vals  = [0, 1, 2]
l2_vals = [0, 2, 6]

# ─────────────────────────────────────────────
# Generate one figure per sector
# ─────────────────────────────────────────────

for (c2_target, Z_even, sector_title, fname) in sectors

    plts = []

    for (col, (L, l2_target)) in enumerate(zip(L_vals, l2_vals))

        p = plot(;
            xlabel     = L"\mathbf{R^{-\omega}}",
            ylabel     = col == 1 ? L"\mathbf{\Delta}" : "",
            title      = L"L = %$L",
            titlefont  = font(11),
            xlims      = (xmin, xmax),
            ylims      = (-0.2, 5.0),
            xticks     = (0.0:0.2:0.6, ["0.0", "0.2", "0.4", "0.6"]),
            yticks     = 0:1:5,
            legend     = false,
            grid       = true,
            gridalpha  = 0.2,
            framestyle = :box,
            size       = (350, 400),
        )

        add_reference_lines!(p, L, c2_target, Z_even)

        # Data points
        sector_data = filter(all_data) do d
            abs(d.l2 - l2_target) < 0.2 + 0.3 * L &&
            iseven(d.Z) == Z_even &&
            abs(d.c2 - c2_target) < 0.5
        end
        println("$(fname), L=$L — nm values: ", sort(unique([d.nm for d in sector_data])))

        for d in sector_data
            scatter!(p, [x_of_nm(d.nm)], [d.Δ];
                markercolor       = nm_to_color[d.nm],
                markerstrokewidth = 0.3,
                markersize        = 5,
                label             = "",
            )
        end

        push!(plts, p)
    end

    # ── Legend panel ─────────────────────────────────────────

    legend_plot = plot(;
        framestyle      = :box,
        axis            = false,
        grid            = false,
        legend          = :topleft,
        legendfontsize  = 9,
        legendrowheight = 20,
        margin          = 0mm,
    )
    for nm in nms
        scatter!(legend_plot, [NaN], [NaN];
            markercolor = nm_to_color[nm],
            markershape = :rect,
            markersize  = 10,
            label       = "nm = $nm",
        )
        scatter!(legend_plot, [NaN], [NaN];
            markercolor       = :white,
            markerstrokewidth = 0,
            markersize        = 0,
            label             = " ",
        )
    end
    scatter!(legend_plot, [NaN], [NaN];
        markercolor = :gray60, markerstrokecolor = :black,
        markershape = :circle, markersize = 6, label = "ED",
    )
    plot!(legend_plot, [NaN, NaN], [NaN, NaN];
        linecolor = :tomato,    linestyle = :solid, linewidth = 1.5,
        label = "descendant (exact)",
    )
    plot!(legend_plot, [NaN, NaN], [NaN, NaN];
        linecolor = :tomato,    linestyle = :dash,  linewidth = 0.5,
        label = "descendant (CB)",
    )
    plot!(legend_plot, [NaN, NaN], [NaN, NaN];
        linecolor = :steelblue, linestyle = :solid, linewidth = 1.5,
        label = "primary (exact)",
    )
    plot!(legend_plot, [NaN, NaN], [NaN, NaN];
        linecolor = :steelblue, linestyle = :dash,  linewidth = 0.5,
        label = "primary (CB)",
    )

    # ── Combine and save ─────────────────────────────────────

    fig = plot(
        plts[1], plts[2], plts[3], legend_plot;
        layout             = @layout([a{0.27w} b{0.27w} c{0.27w} d{0.22w}]),
        size               = (1300, 460),
        plot_title         = sector_title,
        plot_titlefontsize = 18,
        plot_titlevspan    = 0.09,
        margin             = 3mm,
        left_margin        = 8mm,
        right_margin       = 4mm,
        top_margin         = 0mm,
        bottom_margin      = 5mm,
    )

    savefig(fig, joinpath(@__DIR__, fname * ".pdf"))
    savefig(fig, joinpath(@__DIR__, fname * ".png"))
    println("Saved: $fname")
end

println("\nDone — 6 figures generated.")

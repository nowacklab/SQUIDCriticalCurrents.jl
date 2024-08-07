using GLMakie
# Expects definitions from optimize.jl
include("optimize.jl")

function plot_modulation!(ax, parameters; num_phases = 100)
    modulation = critical_current_modulation(parameters; num_phases)
#modulation = upper_critical_current_modulation(parameters; num_phases)

    θas = modulation[:, 1]
    ics = modulation[:, 2]

    positive_ics = 0 .<= ics
    θa_min, θa_max = extrema(θas[positive_ics])
    θa_span = θa_max - θa_min

    locus_color = :gray
    for θa_offset in (-1:1) .* 2π
        lines!(ax, θas .+ θa_offset, ics; color = locus_color)
    end

    for θa_offset in [] .* 2π
        branch_θas = θas[positive_branch]
        perm = sortperm(branch_θas)
        scatterlines!(ax, branch_θas[perm] .+ θa_offset,
                ics[positive_branch][perm]; color = :black)
    end

    return modulation
end

function plot_modulation(parameters; num_phases = 100)
    fig = Figure(;
        title = "$(parameters)",
    )
    ax = Axis(fig[1,1])
    ax.xlabel = rich("Magnetic Flux θ", subscript("a"), " = 2πΦ / Φ", subscript("0"))
    ax.ylabel = rich("Critical current i", subscript("c"),
            " = 2I", subscript("c"), " / (I", subscript("c,1"), " + I", subscript("c,2"),
            ")")
    modulation = plot_modulation!(ax, parameters; num_phases)
    return (fig, ax, modulation)
end

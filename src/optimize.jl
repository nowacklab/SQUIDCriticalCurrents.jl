import Contour
import Roots
using Statistics

# This is a slow implementation for now. Lots of low-hanging fruit.

i1(φ) = sin(φ)
di1(φ) = cos(φ)
i2(φ) = sin(φ)
di2(φ) = cos(φ)

total_current(φ1, φ2; di_0, α1, α2) = (1 + di_0)*i1(φ1) + (1 - di_0)*i2(φ2)
df1(φ1, φ2; di_0, α1, α2) = (1 + di_0)*di1(φ1) + ((1 - di_0)*di2(φ2) / (1 + α2*di2(φ2)))*(1 + α1*di1(φ1))
dfdλ(φ1, φ2; θa, di_0, α1, α2) = φ1 - φ2 + θa + α1*i1(φ1) - α2*i2(φ2)
θa(φ1, φ2; di_0, α1, α2) = φ2 - φ1 - α1*i1(φ1) + α2*i2(φ2)

# di_0 = (ic1 - ic2) / (ic1 + ic2)
parameters = (; di_0 = (1 - 0.8) / (1 + 0.8), α1 = 0.4, α2 = 0.8)

function phase_pairs(parameters; num_phases = 100)
    for iteration in 1:5
        φ1s = φ2s = range(-π, π, length = num_phases)
        dfs = [df1(φ1, φ2; parameters...) for φ1 in φ1s, φ2 in φ2s]
        c = Contour.contour(φ1s, φ2s, dfs, 0.0)
        means = []
        for line in Contour.lines(c)
            xs, ys = Contour.coordinates(line)
            push!(means, mean(xs)^2 + mean(ys)^2)
        end
        cm, i = findmin(means)
        if cm < 1e-2
            xs, ys = Contour.coordinates(Contour.lines(c)[i])
            return [xs ys]
        end
        # Did not sample enough phases to resolve the central contour
        # Try more.
        num_phases *= 4
    end
    throw(ArgumentError("Could not resolve central contour in φ1-φ2 space. Is α1 or α2 much larger than one?"))
end

function interpolate_linear(xs, ys, x)
    if !(length(xs) <= length(ys))
        throw(ArgumentError("The length of xs cannot exceed the length of ys."))
    elseif !(length(xs) >= 2)
        throw(ArgumentError("Need at least two xs to interpolate, but got $(length(xs))"))
    elseif !(xs[1] < xs[end])
        throw(ArgumentError("Interpolating xs must be increasing. Try reverse(xs), reverse(ys)."))
    end
    i = findlast(xs .<= x)
    if isnothing(i) || length(xs) <= i
        if isapprox(x, xs[end])
            return ys[lastindex(xs)]
        elseif isapprox(x, xs[1])
            return ys[firstindex(xs)]
        end
        xmin, xmax = extrema(xs)
        throw(DomainError(x,
            "Cannot interpolate_linear outside the domain [$(xmin), $(xmax)]."))
    end
    x1, x2 = xs[i], xs[i + 1]
    y1, y2 = ys[i], ys[i + 1]
    t = (x - x1) / (x2 - x1)
    return (1 - t)*y1 + t*y2
end

function interpolate_linear_clamped(xs, ys, x)
    return interpolate_linear(xs, ys, clamp(x, xs[1], xs[end]))
end

function interpolate_linear_period(xs, ys, x)
    xmin, xmax = xs[1], xs[end]
    return interpolate_linear(xs, ys, xmin + mod(x - xmin, xmax - xmin))
end

function critical_current_modulation(parameters; num_phases = 100)
    pairs = phase_pairs(parameters; num_phases)
    modulation = [f(pair...; parameters...)
        for pair in eachrow(pairs), f in [θa, total_current]]
    return modulation
end

function upper_critical_current_modulation(modulation)
    "Split critical current locus across line connecting flux-extrema (cusps)"
    θas = modulation[:, 1]
    ics = modulation[:, 2]
    θamin, imin = findmin(θas)
    cuspline(θa) = interpolate_linear([θamin, -θamin], [ics[imin], -ics[imin]], θa)
    upper = [cuspline(θa) <= ic for (θa, ic) in eachrow(modulation)]
    return sortslices(modulation[upper, :], dims = 1)
end

function upper_modulation_endpoints(modulation)
    θas = modulation[:, 1]
    ics = modulation[:, 2]
    θamin = θas[1]
    θamax = -θas[1] - 2π
    domain = (θamin, θamax)
    mrange = θamin .<= θas .<= θamax
    mxs, mys = θas[mrange], ics[mrange]
    fm(θa) = interpolate_linear_clamped(mxs, mys, θa)
    prange = θamin + 2π .<= θas .<= -θamin
    pxs, pys = θas[prange], ics[prange]
    fp(θa) = interpolate_linear_clamped(pxs, pys, θa + 2π)
    θa0 = Roots.find_zero(θa -> fp(θa) - fm(θa), domain)
    return (
        (θa0, fm(θa0)),
        (θa0 + 2π, fp(θa0)),
    )
end

function modulation_endpoints(modulation)
    icp = 0 .<= modulation[:,2]
    θas = modulation[icp, 1]
    ics = modulation[icp, 2]
    θamin, imin = findmin(θas)
    θamax, imax = findmax(θas)

    # Differentiate branch by line between cusps?
    cuspline(θa) = interpolate_linear([θamin, -θamin], [ics[imin], -ics[imin]], θa)
    let xs = range(θamin, -θamin, length = 100)
        lines!(ax, xs, cuspline.(xs))
    end

    θaminrange = θamin .<= θas .<= -2π + θamax
    # TODO: Check that this Δθa sign heuristic always works.
    # Nope.
    finalsign = sign(θas[end] - θas[end-1])
    dθarange = finalsign .== sign.([[θas[1] - θas[end]]; θas[2:end] .- θas[1:end-1]])
    θamaxrange = 2π + θamin .<= θas .<= θamax
    minrange = θaminrange .&& dθarange
    mxs, mys = θas[minrange], ics[minrange]
    mperm = sortperm(mxs)
    mxs, mys = mxs[mperm], mys[mperm]
    maxrange = θamaxrange .&& dθarange
    pxs, pys = θas[maxrange], ics[maxrange]
    pperm = sortperm(pxs)
    pxs, pys = pxs[pperm], pys[pperm]
    domain = (max(mxs[1], pxs[1] - 2π), min(mxs[end], pxs[end] - 2π))
    fm(θa) = interpolate_linear(mxs, mys, θa)
    fp(θa) = interpolate_linear(pxs, pys, θa + 2π)
    θa0 = Roots.find_zero(θa -> fp(θa) - fm(θa), domain)
    return (
        (θa0, fm(θa0)),
        (θa0 + 2π, fp(θa0)),
    )
end

function squid_has_vanishing_critical_current(parameters)
    (;di_0, α1, α2) = parameters
    return  di_0 == 0 && α1 == 0 && α2 == 0
end

function upper_modulation(parameters; num_phases = 100)
    if squid_has_vanishing_critical_current(parameters)
        θas = range(-π, π, length = 4*num_phases)
        ics = [2*abs(cos(θa / 2)) for θa in θas]
        return (θas, ics)
    end
    modulation = critical_current_modulation(parameters; num_phases)
    θas = modulation[:, 1]
    ics = modulation[:, 2]

    positive_ics = 0 .<= ics
    pθas = θas[positive_ics]
    θa_min, θa_max = extrema(pθas)
    θa_span = θa_max - θa_min
    if θa_span <= 2π
        pics = ics[positive_ics]
        perm = sortperm(pθas)
        return (pθas[perm], pics[perm])
    end
    up_modulation = upper_critical_current_modulation(modulation)
    up_θas = up_modulation[:, 1]
    up_ics = up_modulation[:, 2]
    (θa_min, ic_min), (θa_max, ic_max) = upper_modulation_endpoints(up_modulation)
    inperiod = θa_min .< up_θas .< θa_max
    return (
        [[θa_min]; up_θas[inperiod]; [θa_max]],
        [[ic_min]; up_ics[inperiod]; [ic_max]],
    )
end

function positive_modulation(parameters; num_phases = 100)
    if squid_has_vanishing_critical_current(parameters)
        θas = range(-π, π, length = 4*num_phases)
        ics = [2*abs(cos(θa / 2)) for θa in θas]
        return (θas, ics)
    end
    modulation = critical_current_modulation(parameters; num_phases)
    θas = modulation[:, 1]
    ics = modulation[:, 2]

    positive_ics = 0 .<= ics
    pθas = θas[positive_ics]
    θa_min, θa_max = extrema(pθas)
    θa_span = θa_max - θa_min
    if θa_span <= 2π
        pics = ics[positive_ics]
        perm = sortperm(pθas)
        return (pθas[perm], pics[perm])
    end
    (θa_min, ic_min), (θa_max, ic_max) = modulation_endpoints(modulation)
    finalsign = sign(pθas[end] - pθas[end-1])
    dθarange = finalsign .== sign.([[θas[1] - θas[end]]; θas[2:end] .- θas[1:end-1]])
    positive_branch = positive_ics .&& dθarange .&& θa_min .<= modulation[:, 1] .<= θa_max
    branch_θas = modulation[:, 1][positive_branch]
    branch_ics = modulation[:, 2][positive_branch]
    perm = sortperm(branch_θas)
    return (
        [[θa_min]; branch_θas[perm]; [θa_max]],
        [[ic_min]; branch_ics[perm]; [ic_max]],
    )
end

nothing;


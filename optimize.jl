import Contour
import Roots
using Statistics

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
    φ1s = φ2s = range(-π, π, length = num_phases)
    dfs = [df1(φ1, φ2; parameters...) for φ1 in φ1s, φ2 in φ2s]
    c = Contour.contour(φ1s, φ2s, dfs, 0.0)
    coordinates = [Contour.coordinates(line) for line in Contour.lines(c)]
    allxs = empty(φ1s)
    allys = empty(φ2s)
    for coords in coordinates
        xs, ys = coords
        append!(allxs, xs)
        append!(allys, ys)
    end
    return [allxs allys]
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
        if x == xs[end]
            return ys[lastindex(xs)]
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

function interpolate_linear_period(xs, ys, x)
    xmin, xmax = xs[1], xs[end]
    return interpolate_linear(xs, ys, xmin + mod(x - xmin, xmax - xmin))
end

function critical_current_modulation(parameters; num_phases = 100)
    pairs = phase_pairs(parameters; num_phases)
    modulation = [f(pair...; parameters...)
        for pair in eachrow(pairs), f in [θa, total_current]]
end

function modulation_endpoints(modulation)
    # TODO: Sometimes the curve in ics is not single-valued (at the cusp).
    # Informally, the intersection on the positive branch is above the cusp.
    # This is breaking the endpoint detection.
    icp = 0 .<= modulation[:,2]
    θas = modulation[icp, 1]
    ics = modulation[icp, 2]
    θamin, θamax = extrema(θas)
    minrange = θamin .<= θas .<= -2π + θamax
    maxrange = 2π + θamin .<= θas .<= θamax
    mxs, mys = ics[minrange], θas[minrange]
    fm = x -> interpolate_linear(mxs, mys, x)
    pxs, pys = reverse(ics[maxrange]), reverse(θas[maxrange])
    fp = x -> interpolate_linear(pxs, pys, x)
    domain = (max(mxs[1], pxs[1]), min(mxs[end], pxs[end]))
    ic0 = Roots.find_zero(x -> fp(x) - fm(x) - 2π, domain)
    return ((fm(ic0), ic0), (fp(ic0), ic0))
end

function squid_has_vanishing_critical_current(parameters)
    (;di_0, α1, α2) = parameters
    return  di_0 == 0 && α1 == 0 && α2 == 0
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
    positive_branch = positive_ics .&& θa_min .<= modulation[:, 1] .<= θa_max
    branch_θas = modulation[:, 1][positive_branch]
    branch_ics = modulation[:, 2][positive_branch]
    perm = sortperm(branch_θas)
    return (
        [[θa_min]; branch_θas[perm]; [θa_max]],
        [[ic_min]; branch_ics[perm]; [ic_max]],
    )
end

nothing;


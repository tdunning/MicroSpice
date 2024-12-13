using MicroSpice

include("circuits.jl")
include("components.jl")
include("optimize.jl")
include("plots.jl")

# fundamental for 10m WSPR transmission
f = 28.126100e6

# harmonic suppression for 55% duty cycle square wave
r0 = [-19,-14.5,-24,-27]

function q()
    sim = MicroSpice.solve(nl6)
    params -> quality(sim, params, f:f:5f, r0)
end

components = [inductors(true), caps(), caps(), caps()]
start = rand.(eachindex.(components))

h,best = metaClimb(q, start, components, gens=100, rate=1, λ=0.1, pop=100, μ=0.9, info=10)

# check out the 5 best performers at the end
best[1:5]

# plot best performers over optimization process
statusPlot(nl6, getindex.(components, h[2][1].ix), f:f:5f, r0)
statusPlot!(nl6, getindex.(components, h[10][1].ix), f:f:5f, r0)
statusPlot!(nl6, getindex.(components, h[20][1].ix), f:f:5f, r0)
statusPlot!(nl6, getindex.(components, h[50][1].ix), f:f:5f, r0)
gui()

# plot the size of the mutation steps
using LinearAlgebra
plot(mean(norm.(getfield.(getindex.(h, (1:10)'), :ix) .- getfield.(getindex.(h, (1:10)'), :old)), dims=2))

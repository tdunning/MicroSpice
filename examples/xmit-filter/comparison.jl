# This sc
using CSV, DataFrames
using Interpolations

ref = CSV.read("examples/xmit-filter/images/reference.csv", DataFrame)

c = linear_interpolation(ref[:,:f], ref[!,:response])

statusPlot(nl6, getindex.([inductors(false), caps(), caps(), caps()], h[100][1].ix), f:f:5f, r0)
plot!(ref[!,:f]/1e6, ref[!,:response] .+ 20*log10(3/2), label="Designed")
plot!(size=(800,600), xlim=(0,200), legend=true)
scatter!(f->c(f*1e6)+20*log10(3/2), (f:f:5f) ./ 1e6, label=false)

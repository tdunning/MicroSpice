using Plots

"""
Plots the response curve for a model (a function of frequency
returning dB response) with highlights on the frequencies of interest
and table showing how the figure of merit for a filter is
computed. The frequencies given are a vector of harmonics starting
with the fundamental. The penalties correspond to the degree of
inherent suppression of the harmonics for the 2nd through nth
harmonics.
"""
function statusPlot(nl, parameters, frequencies, penalties)
    plot(legend=false, yrange=(-60, 5),
         xlab="Frequency (MHz)\n", ylab="\nResponse (dB)",
         title="\nFilter analysis", size=(900, 700))
    statusPlot!(nl, parameters, frequencies, penalties)
end

function annotatePlot!(model, frequencies, penalties)
    y = -20
    dy = 1.5
    x0 = 200
    i = 1
    annotate!(x0, y+dy, ("Response", :right, 10))
    annotate!(x0 + 30, y+dy, ("Net", :right, 10))
    response = model(frequencies[1])
    r1 = response
    annotate!(x0, y, ("$(round(response, digits=1))", :right, 8))
    annotate!(x0 + 30, y, ("--", :right, 8))
    y -= dy
    for f in frequencies[2:end]
        response = model(f)
        annotate!(x0, y, ("$(round(response, digits=1))", :right, 8))
        annotate!(x0 + 30, y, ("$(round(response - r1, digits=1))", :right, 8))
        annotate!(x0 + 60, y, ("$(round(response - r1 + penalties[i], digits=1))", :right, 8))
        y -= dy
        i += 1
    end
end

"""
Adds an additional plot
"""
function statusPlot!(nl, parameters, frequencies, penalties)
    expRange(low, high, n) = exp.(LinRange(log(low), log(high), n))

    fx = MicroSpice.solve(nl)
    model = f->20*log10(abs(only(fx(f, [2,0], parameters))))
    
    top = model(frequencies[1]) + 3
    plot!(f -> model(f * 1e6), expRange(10,300,100), label="Evolved")
    scatter!(f -> model(f * 1e6), frequencies ./ 1e6, markercolor="green", label=false)
    scatter!(frequencies[2:end]./1e6, model.(frequencies[2:end]) .+ penalties, markercolor="red", label="Output level")
    return model
end

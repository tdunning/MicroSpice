using MicroSpice

nl1 = MicroSpice.Netlist(raw"""
R1 in N001 50
L1 N001 out $L1
C1 N001 gnd $C1
C2 out N002 $C2
R2 out gnd 50
L2 N002 gnd $L2
R3 N002 gnd $R
""", [:L1, :C1, :C2, :L2, :R])

nl2 = MicroSpice.Netlist(raw"""
R1 in N001 50
L1 N001 out $L1
C1 N001 gnd $C1
C2 out gnd $C2
R2 out gnd 50
""", [:L1, :C1, :C2])
c2 = [300e-9, 150e-12, 150e-12]

nl3 = MicroSpice.Netlist(raw"""
R1 in out 50
C1 out N001 225p
R2 out gnd 50
L2 N001 gnd 20n
R3 N001 gnd 30
""", [])

"""
Chebyshev 5th order with 5db passband ripple
"""
nl4 = MicroSpice.Netlist(raw"""
R1 in N001 50
L1 N001 N002 $L1
C1 N002 0 $C2
R2 out 0 50
C2 N001 0 $C1
L2 N002 out $L2
C3 out 0 $C3
""", [:L1, :L2, :C1, :C2, :C3])
c4 = [162.1753e-9, 162.1753e-9, 550.6910e-12, 700.1618e-12, 550.6910e-12]

"""
Chebyshev 3rd order with 5dB passband ripple
"""
nl5 = MicroSpice.Netlist(raw"""
R1 in N001 50
L1 N001 out $L1
V1 in 0 SINE() AC 2
C2 out 0 $C2
R2 out 0 50
C1 N001 0 $C1
""", [:L1, :C1, :C2])
c5 = [151.9541e-9, 529.6833e-12, 529.6833e-12]

"""Like LinRange, but returns values evenly spaced on a log-scale"""
expRange(low, high, n) = exp.(LinRange(log(low), log(high), n))
                              
"""
Plots the response curve for a model (a function of frequency
returning dB response) with highlights on the frequencies of interest
and table showing how the figure of merit for a filter is
computed. The frequencies given are a vector of harmonics starting
with the fundamental. The penalties correspond to the degree of
inherent suppression of the harmonics for the 2nd through nth
harmonics.
"""
function statusPlot(model, frequencies, penalties)
    modelx = f -> model(f*1e6)
    
    top = model(frequencies[1]) + 3
    plot(modelx, expRange(10,300,100), legend=false, yrange=(top-35, top),
         xlab="Frequency (MHz)\n", ylab="\nResponse (dB)",
         title="\nFilter analysis")
    scatter!(modelx, frequencies ./ 1e6)
    scatter!(frequencies[2:end]./1e6, model.(frequencies[2:end]) .- penalties)
    y = top-20
    dy = 1.7
    x0 = 250
    i = 1
    annotate!(x0, y+dy, ("Response", :right, 10))
    annotate!(x0 + 30, y+dy, ("Net", :right, 10))
    response = model(frequencies[1])
    annotate!(x0, y, ("$(round(response, digits=1))", :right, 8))
    annotate!(x0 + 30, y, ("--", :right, 8))
    y -= dy
    for f in frequencies[2:end]
        response = model(f)
        annotate!(x0, y, ("$(round(response, digits=1))", :right, 8))
        annotate!(x0 + 30, y, ("$(round(response-penalties[i], digits=1))", :right, 8))
        y -= dy
        i += 1
    end
end

function fitModel(nl, x0, frequencies, penalties)
    MicroSpice.solve(nl, [:in, :gnd], [:out], parameters)
    model = f->20*log10(abs(only(fx(f, [2,0]))))
    q = maximum(merit(model, frequencies, penalties))
end

function merit(model, frequencies, penalties)
    response = model.(frequencies)
    net = (response[1] .- response[2:end]) .- penalties
end
    
function resample(nl, quality, samples=20)
    parameters -> begin
        p = (0.95 .+ 0.1 * rand(samples, length(parameters))) .* parameters'
        minimum([quality(nl, c) for c in eachrow(p)])
    end
end

function attenuation(frequencies, penalties)
    offset = 10 .^ (-penalties/20)
    (nl, parameters) -> begin
        s = MicroSpice.solve(nl, [:in, :gnd], [:out], parameters)
        r = (abs ∘ only ∘ s).(frequencies, Ref([1,0]))
        r[1].^2 / sum(offset .* r[2:end])
    end
end

f = 28.126100e6
r0 = [-19,-14.5,-24,-27]
base = [220e-9, 220e-12, 220e-12, 10e-9, 30]
base2 = [220e-9, 220e-12, 220e-12, 40e-9, 100]
#frequencies = f * (1:5)
decibel(x) = 20 * log10(abs(x))

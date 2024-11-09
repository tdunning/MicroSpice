module MicroSpice

using LinearAlgebra

# A collection of simple functions useful in electronic design

"Inductance in nH of a coil with n turns d mm in diameter and l mm long"
coil(n, d, l) = (l * n)^2 / (0.44 * d + l)

"Complex reactance of an inductor at a given frequency"
xl(l, f) = 2π*f*l * 1.0im

"Complex reactance of a capacitor at a given frequency"
xc(c, f) = 1/(2π*f*c * 1.0im)

"""
A Link is a composite of pure resistance, inductance and capacitance
arranged in parallel. 

To simplify combining multiple links, resistance and inductance are 
retained in inverse form (i.e. 1/R, 1/L).
"""
struct Link
    r
    c
    l
end
Base.isapprox(x::Link, y::Link) = x.r ≈ y.r && x.c ≈ y.c && x.l ≈ y.l

Cap(c) = Link(0, c, 0)
Ind(L) = Link(0, 0, 1.0/L)
Res(r) = Link(1/r, 0, 0)

Base.:+(x::Link, y::Link) = Link(x.r + y.r, x.c + y.c, x.l + y.l)

impedance(x::Link, f) = 1/(1/x.r + 1/xc(x.r, f) + 1/xl(x.l, f))

"""
A netlist is defined as a set of nodes which are connected by
ideal resistors, inductors and capacitors. Such a netlist can
be resolved into a ComplexCircuit at any desired frequency by
replacing the inductors and capacitors by their equivalent 
complex impedances at that frequency.
"""
struct Netlist
    n::Int
    links::Matrix{Link}
    names::Dict{String,Int}
end


"""
A circuit is defined as a set of nodes connected by admittances. Each
node may be set to a particular voltage, or can have a current injected
or it can be free-floating. This is modeled as a linear system of equations 
in terms of the `n` nodal voltages of which `m` are externally constrained, 
and `k` injected currents. The first `n` of these equations define 
the relationship between node voltages and injected currents and are
of rank `n-1`. The next `m` equations define externally set voltages 
followed by `k` equations for injected currents. We assume that imposed
voltages and injected currents are on different nodes so this gives
a rank of `n + m + k - 1` so far. There are an additional `n-m-k`
equations that set all remaining injected currents to zero and a
final constraint that sets the sum of all injected currents to zero.
"""
struct ComplexCircuit
    n::Int
    links::Matrix
end


raw"""
Constructs a Netlist from a Spice circuit. As an example, here 
is a simple RLC circuit with `in` as input and `out` as output
```
L1 in  out 100nH
R1 out gnd  50Ω
C1 out gnd 100nF
```
Currently the only supported components are resistors (R),
inductors (L) and capacitors (C). Component values can be
specified using standard SI prefixes such as `k` (kilo = 1e3), 
`M` (Mega=1e6), `G`(giga=1e9), `m` (milli=1e-3), `μ` (micro=1e-6), 
`n` (nano=1e-9), or `p` (pico=1e-12). In addition, the non-standard 
forms `u` (micro=1e-6) and `meg` (Mega=1e6) are used for Spice
compatibility.

Units including H (Henry), F (Farad) and Ω (ohm) are allowed, but 
are optional, but only the correct unit type is allowed.f

A Netlist can be converted to a ComplexCircuit at a specific frequency
and then solved with specified input voltages and currents to find
the voltages for all nodes in the circuit. Taking the example above, 
we can find the response for a few frequencies using this 
code
```jldoctest; filter = r"(\d*)\.(\d{9})\d+" => s"\1.\2***"
nl = MicroSpice.Netlist("L1 in  out 100n\nR1 out gnd 50\nC1 out gnd 100n\n")
s = MicroSpice.solve(nl, [:in, :gnd], [:out])
decibel(x) = 20 * log10(abs(x))
[decibel(only(s(f, [1, 0]))) for f in [1.4e6, 1.5e6, 1.62e6, 1.8e6]]
# output
4-element Vector{Float64}:
 12.883077832402897
 18.914298710446417
 27.65586909592151
 11.05634871505561
```
"""
Netlist(s::AbstractString) = Netlist(IOBuffer(s))
function Netlist(input)
    c = []
    nodes = Dict()
    links = Dict{Tuple,Link}()
    for x in readlines(input)
        kind = uppercase(x[1])
        if !(kind in "RLC")
            continue
        end
        name, from, to, v = split(x)
        value = decode(kind, v)
        i = nodes[from] = get(nodes, from, length(nodes)+1)
        j = nodes[to] = get(nodes, to, length(nodes)+1)
        i, j = sort([i, j])
        
        if kind == 'C'
            z = Cap(value)
        elseif kind == 'L'
            z = Ind(value)
        elseif kind == 'R'
            z = Res(value)
        else
            continue
        end
        links[(i, j)] = get(links, (i, j), Link(0, 0, 0)) + z
    end
    n = length(nodes)
    m = fill(Link(0,0,0), (n,n))
    for (coords,v) in links
        m[coords...] = v
        m[reverse(coords)...] = v
    end
    Netlist(n, m, nodes)
end

function solve(nl::Netlist, inputs::Vector, outputs::Vector)
    vIn = [nl.names[string(k)] for k in inputs]
    vOut = [nl.names[string(k)] for k in outputs]
    (frequency, inputs) -> begin
        c = ComplexCircuit(nl, frequency)
        result = solve(c, zip(vIn,inputs), [])
        return getindex.(Ref(result), vOut)
    end
end

function inputImpedance(nl::Netlist, in, ref)
    input = nl.names[string(in)]
    gnd = nl.names[string(ref)]
    frequency -> begin
        c = ComplexCircuit(nl, frequency)
        state = solve(c, [gnd => 0], [input => 1])
        return state[input]
    end
end

function transfer(nl::Netlist, in, out, ref)
    input = nl.names[string(in)]
    output = nl.names[string(out)]
    gnd = nl.names[string(ref)]
    frequency -> begin
        c = ComplexCircuit(nl, frequency)
        state = solve(c, [input => 1, gnd => 0], [])
        return state[output]
    end
end


function decode(kind, x::AbstractString)
    scaleFactor = Dict("k"=>1e3, "M"=>1e6, "meg"=>1e6, "mega"=>1e6, "G"=>1e9, 
                       "m"=>1e-3, "u"=>1e-6, "μ"=>1e-6, "n"=>1e-9, "p"=>1e-12,
                       ""=>1.0, nothing=>1.0)
    units = Dict("H" => 'L', "F" => 'C', "Ω" => 'R', "ohm" => 'R')

    m = match(r"(\d+(?:\.\d*)?(?:e\d+)?)\s*((?:mega|meg)|[mkMGμunp])?\s*((?:ohm)|[HFΩ])?", x)
    if m === nothing
        return 0
    else
        raw, scale, unit = m
        if unit !== nothing && unit != "" && kind ≠ units[unit]
            throw(ErrorException("Bad unit ($unit) for component kind $kind"))
        end
        return parse(Float64, raw) * scaleFactor[scale]
    end
end



size(c::ComplexCircuit) = c.n

function ComplexCircuit(n::Netlist, f)
    σ = admittance.(n.links, f)
    i = diagind(σ)
    σ[i] .= 0
    σ[i] = -sum(σ, dims=1)
    ComplexCircuit(length(n.names), σ)
end

admittance(x::Link, f) = let ω = 2π * f * -1im
    x.r + ω * x.c + x.l / ω
end

"""
Compute the voltages and injected currents for each node.

We write node voltages as `vᵢ`, injected currents as `Iᵢ`,
and node to node complex admittance as `σᵢⱼ`. To 
express the node-to-node currents in matrix form, we set
the diagonal of `σ` to the negative sum of the row.

Internally, we define linear equations in terms of the `n` node
voltages of which `k` are set and `n` injected node currents all
but `k` of which are zero or known. We solve for all `n` voltages
and `n` injected currents.

External voltage and injected current is given a list of 
pairs of node name and voltage or current.

The return value is a vector of node voltages.
"""
function solve(c::ComplexCircuit, externalVoltages, injectedCurrents) 
    n = size(c)
    
    # relationship of voltages and injected currents
    # \sum_{j\ne i} (v_j - v_i) \sigma_{ij} - I_i = 0
    # but \sigma_{ii} = -\sum_{j \ne i} \sigma_{ij} so
    # \sum_{j=1\ldots n} \sigma_{ij} v_i - I_i = 0
    eq1 = [c.links -I]
    rhs1 = zeros(n)

    # impose the voltages (if any).
    eq2 = [I zeros(n, n)][first.(externalVoltages), :]
    rhs2 = getindex.(externalVoltages, 2)

    # set the net currents to zero, 
    eq3 = [zeros(n, n) I]
    rhs3 = zeros(n)

    # except where we inject currents
    i = first.(injectedCurrents)
    rhs3[i] = -getindex.(injectedCurrents, 2)

    # and remove equations where we set the voltage and
    # thus don't know the net current
    i = setdiff(1:n, first.(externalVoltages))
    eq3 = eq3[i,:]
    rhs3 = rhs3[i]

    # and finally, all injected (and other) currents must balance
    eq4 = [zeros(n)' ones(n)']
    rhs4 = zeros(1)

    eq = [eq1; eq2; eq3; eq4]
    rhs = [rhs1; rhs2; rhs3; rhs4]
    return eq \ rhs
end


# function r2r(n)
#     f = ComplexCircuit(n)
#     for i in 2:(n-1)
#         addLink(f, i, 1, 2)
#     end
#     addLink(f, n, 1, 1)
#     for i in 3:n
#         addLink(f, i, i-1, 1)
#     end
#     f
# end
# 
# 
# function rc(f)
#     c = ComplexCircuit(3, Matrix{Complex{Num}}(zeros(3,3)))
#     addLink(c, 1, 2, xc(1e-9, f))
#     addLink(c, 2, 3, 1e3)
#     c
# end
# 
# using Plots, Symbolics, LinearAlgebra
# 
# @variables f
# 
# c = rc(f)
# gain(fx) = solve(c, [(1,0), (3, 1)], [], Dict(f => fx))[2]
# plot(fx -> 20*log10(abs(gain(fx))), [1,2,5,10,20,50,100,200,500,1000,2000,5000,10000] .* 1e3, xscale=:log10)
# 

end

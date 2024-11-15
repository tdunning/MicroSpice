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
An InternalExpression is a sum of a numerical value and references
into a value vector that will be provided at simulation time.
"""
struct InternalExpression
    invert::Bool
    base::Float64
    refs::Vector{Int}

    function InternalExpression(inv::Bool, x::Number, names::Vector)
        new(inv, x, names)
    end
end

InternalExpression(kind::AbstractChar, x::InternalExpression) = x

function InternalExpression(kind::AbstractChar, x::Number, names::Vector)
    InternalExpression(kind, x)
end

function InternalExpression(kind::AbstractChar, x::Number)
    InternalExpression(kind in "RL", x, [])
end

function InternalExpression(kind::AbstractChar, x::AbstractString, names::Vector)
    n = findfirst(z->String(z)==x, names)
    if isnothing(n)
        throw(ErrorException("Can't find $x in known names: $names"))
    end
    InternalExpression(kind in "RL", 0, [n])
end

function InternalExpression(kind::AbstractChar, x::Symbol, names::Vector)
    InternalExpression(kind, String(x), names)
end

function Base.:+(x::InternalExpression, y::InternalExpression)
    if x.invert != y.invert
        throw(ErrorException("Incompatible values: $x + $y"))
    end
    InternalExpression(x.invert, x.base + y.base, vcat(x.refs, y.refs))
end

function Base.isapprox(x::InternalExpression, y::InternalExpression; kwargs...)
    isapprox(x.base, y.base; kwargs...) && Set(x.refs) == Set(x.refs)
end

function resolve(x::InternalExpression, values::AbstractVector)
    if length(x.refs) == 0
        # special case here to avoid sum(Vector{Any}[]) propagating bad type info
        return x.base
    else
        if length(values) == 0
            throw(BoundsError("No values provided for reference value"))
        end
        raw = getindex.(Ref(values), x.refs)
        if x.invert
            f = x->1/x
            x.base + sum(f.(raw))
        else
            x.base + sum(raw)
        end
    end
end

"""
A Link is a composite of pure resistance, inductance and capacitance
arranged in parallel. 

To simplify combining multiple links, resistance and inductance are 
retained in inverse form (i.e. 1/R, 1/L).
"""
struct Link
    r::InternalExpression
    c::InternalExpression
    l::InternalExpression
    function Link(rx, cx, lx)
        new(InternalExpression('R', rx), InternalExpression('C', cx), InternalExpression('L', lx))
    end

    function Link(r, c, l, names::Vector)
        Link(InternalExpression(r, names), 
             InternalExpression(l, names),
             InternalExpression(c, names))
    end
end


Base.isapprox(x::Link, y::Link) = x.r ≈ y.r && x.c ≈ y.c && x.l ≈ y.l

# For actual values, we invert now, if necessary
Cap(C::Number) = Link(  0, C,   0)
Ind(L::Number) = Link(  0, 0, 1/L)
Res(R::Number) = Link(1/R, 0,   0)

# For expressions, the inversion is already done or will happen at resolve time
Cap(C::InternalExpression) = Link(0, C, 0)
Ind(L::InternalExpression) = Link(0, 0, L)
Res(R::InternalExpression) = Link(R, 0, 0)

Base.:+(x::Link, y::Link) = Link(x.r + y.r, x.c + y.c, x.l + y.l)

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


"""
A netlist is defined as a set of nodes which are connected by
ideal resistors, inductors and capacitors. Such a netlist can
be resolved into a ComplexCircuit at any desired frequency by
replacing the inductors and capacitors by their equivalent 
complex impedances at that frequency.

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
```jldoctest; filter = r"(\\d*)\\.(\\d{9})\\d+" => s"\\1.\\2***"
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
Netlist(s::AbstractString, names=[]) = Netlist(IOBuffer(s), names)

"""
Reads some input and parses that input as a netlist
"""
function Netlist(input, names=[])
    c = []
    nodes = Dict()
    links = Dict{Tuple,Link}()
    for x in readlines(input)
        kind = uppercase(x[1])
        if !(kind in "RLC")
            continue
        end
        name, from, to, v = split(x)
        value = decode(kind, v, names)
        i = nodes[from] = get(nodes, from, length(nodes)+1)
        j = nodes[to] = get(nodes, to, length(nodes)+1)
        i, j = sort([i, j])
        
        a = get(links, (i, j), Link(0, 0, 0))
        links[(i, j)] = a + value
    end
    n = length(nodes)
    m = fill(Link(0, 0, 0), (n, n))
    for (coords,v) in links
        m[coords...] = v
        m[reverse(coords)...] = v
    end
    Netlist(n, m, nodes)
end

"""
Returns a function that simulates a circuit with specified voltage
inputs and outputs. That function takes as arguments the frequency
and a list of input voltages and it returns a list of output voltages.
"""
function solve(nl::Netlist, inputs::AbstractVector, outputs::AbstractVector, parameters::AbstractVector=[])
    vIn = [nl.names[string(k)] for k in inputs]
    vOut = [nl.names[string(k)] for k in outputs]
    (frequency, inputs) -> begin
        c = ComplexCircuit(nl, frequency, parameters)
        result = solve(c, zip(vIn,inputs), [])
        return getindex.(Ref(result), vOut)
    end
end

""" 
Returns a function that takes a frequency and computes the
impedance of a circuit from the specified input to the specified
reference (ground) node.
"""
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

function decode(kind, x::AbstractString, names::Vector=[])
    scaleFactor = Dict("k"=>1e3, "M"=>1e6, "meg"=>1e6, "mega"=>1e6, "G"=>1e9, 
                       "m"=>1e-3, "u"=>1e-6, "μ"=>1e-6, "n"=>1e-9, "p"=>1e-12,
                       ""=>1.0, nothing=>1.0)
    units = Dict("H" => 'L', "F" => 'C', "Ω" => 'R', "ohm" => 'R')

    if x[1] == '$'
        value = InternalExpression(kind, x[2:end], names)
    else
        m = match(r"(\d+(?:\.\d*)?(?:e\d+)?)\s*((?:mega|meg)|[mkMGμunp])?\s*((?:ohm)|[HFΩ])?", x)
        if m === nothing
            value = 0
        else
            raw, scale, unit = m
            if unit !== nothing && unit != "" && kind ≠ units[unit]
                throw(ErrorException("Bad unit ($unit) for component kind $kind"))
            end
            value = parse(Float64, raw) * scaleFactor[scale]
        end
    end
    if kind == 'C'
        return Cap(value)
    elseif kind == 'L'
        return Ind(value)
    elseif kind == 'R'
        return Res(value)
    else
        throw(ErrorException("Bad kind of component: $kind"))
    end
end


size(c::ComplexCircuit) = c.n

function ComplexCircuit(n::Netlist, f, parameters=[])
    let ω = 2π * f * -1im,
        f = link -> 
            let lx = resolve.([link.r, link.l, link.c], Ref(parameters))
                r, l, c = lx
                r + ω * c + l / ω
            end
        σ = f.(n.links)
    i = diagind(σ)
    σ[i] = -sum(σ, dims=1)
    ComplexCircuit(length(n.names), σ)
end

end

"""
Compute the voltages and injected currents for each node.

The input is a list of ``(i, V_i)`` pairs which represent the
imposed voltages and a list of injected currents as ``(i, I_i)``.
The nodes with imposed voltages should not overlap with the nodes
with injected currents. All nodes without an imposed voltage are
considered free-floating and their voltage level will be computed.
Nodes with no injected current are considered isolated except for
the connections specified in the connections in the circuit.

The result is a vector of node voltages and injected currents
``[V_1 \\ldots V_n, I_1 \\ldots I_n]``

"""
function solve(c::ComplexCircuit, externalVoltages, injectedCurrents) 
    n = size(c)
    
    # relationship of voltages and injected currents
    # \sum_{j\ne i} (v_j - v_i) \sigma_{ij} - I_i = 0
    # but \sigma_{ii} = -\sum_{j \ne i} \sigma_{ij} so
    # \sum_{j=1\ldots n} \sigma_{ij} v_i - I_i = 0

    # we construct an identity matrix the hard way to ensure type matching
    function id(n)
        v = zeros(typeof(c.links[1,1]), n, n)
        v[diagind(v)] .= 1
        return v
    end

    t = typeof(c.links[1,1])
        
    eq1 = [c.links -id(n)]
    rhs1 = zeros(t, n)

    # impose the voltages (if any).
    eq2 = [id(n) zeros(t, n, n)][first.(externalVoltages), :]
    rhs2 = zeros(t, length(externalVoltages), 1)
    rhs2 .= getindex.(externalVoltages, 2)

    # set the net currents to zero, 
    eq3 = [zeros(t, n, n) id(n)]
    rhs3 = zeros(t, n, 1)

    # except where we inject currents
    i = first.(injectedCurrents)
    rhs3[i] = -getindex.(injectedCurrents, 2)

    # and remove equations where we set the voltage and
    # thus don't know the net current
    i = setdiff(1:n, first.(externalVoltages))
    eq3 = eq3[i,:]
    rhs3 = rhs3[i]

    # and finally, all injected (and other) currents must balance
    eq4 = [zeros(t, n)' ones(t, n)']
    rhs4 = zeros(t, 1)

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

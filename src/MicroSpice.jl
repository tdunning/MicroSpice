module MicroSpice

using LinearAlgebra


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

"""
A circuit is defined as a set of nodes connected by admittances. Each
node may be set to a particular voltage, or can have a current injected
or it can be free-floating. This is modeled as a linear system of equations
in terms of the `n` nodal voltages of which `m` are externally constrained,
and `k` injected currents. The first `n` of these equations define
the relationship between node voltages and injected currents and have
of rank `n-1`. The next `m` equations define externally set voltages
followed by `k` equations for injected currents. We assume that imposed
voltages and injected currents are on different nodes so this gives
a rank of `n + m + k - 1` so far. There are an additional `n-m-k`
equations that set all remaining injected currents to zero and a
final constraint that sets the sum of all injected currents to zero.

A `Netlist` is normally constructed by parsing a Spice compatible
netlist that defines the ideal components that connect the nodes in
the circuit. The first step in solving a circuit is to bind values to
any symbolically defined component values in the circuit. The second
step is to define a frequency of interest and the circuit inputs to
get a vector of outputs. Because the underlying matrices are part of
the `Netlist` structure itself, solving a circuit is not thread-safe.
"""
struct Netlist
    n::Int                   # number of nodes in the circuit
    links::Dict{Tuple, Link} # abstract admittance matrix (n x n) in sparse form
    names::Dict{String,Int}  # mapping from node name to node number
    inputs::Vector           # indexes of input nodes
    outputs::Vector          # indexes of output nodes 
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

Solving proceeds in two steps. In the first step, any symbolically
defined component values are bound to specific values to get a
function to solve the bound circuit. In the second step, the circuit
is solved for a specific frequency and vector of input voltages.

```jldoctest; filter = r"(\\d*)\\.(\\d{9})\\d+" => s"\\1.\\2***"
nl = MicroSpice.Netlist("L1 in  out 100n\nR1 out gnd 50\nC1 out gnd 100n\n", [], [:in, :gnd], [:out])
s = MicroSpice.solve(nl)
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
Netlist(s::AbstractString, names, inputs, outputs) = Netlist(IOBuffer(s), names, inputs, outputs)

"""
Reads some input and parses that input as a netlist
"""
function Netlist(input, parameter_names, inputs, outputs)
    c = []
    nodes = Dict()
    links = Dict{Tuple,Link}()
    n = 0
    for x in readlines(input)
        kind = uppercase(x[1])
        if !(kind in "RLC")
            continue
        end
        name, from, to, v = split(x)
        value = decode(kind, v, parameter_names)
        i = nodes[from] = get(nodes, from, length(nodes)+1)
        j = nodes[to] = get(nodes, to, length(nodes)+1)
        i, j = sort([i, j])
        n = max(n, i, j)

        a = get(links, (i, j), Link(0, 0, 0))
        links[(i, j)] = a + value
    end
    inputNodes = getindex.(Ref(nodes), string.(inputs))
    outputNodes = getindex.(Ref(nodes), string.(outputs))

    Netlist(n, links, nodes, inputNodes, outputNodes)
end

"""
Given a `Netlist` return a function that will solve the circuit
for all voltages and injected currents given frequency, inputs
and circuit parameters. 

The returned value is a vector of the n output voltages and n injected
currents.
"""
function raw(nl::Netlist)
    # the system of equations is allocated once, solved many times
    n = nl.n
    t = Complex{Float64}
    eq1 = [zeros(t, n, n) -I]
    rhs1 = zeros(t, n)
    
    # impose input voltages (if any).
    eq2 = [I zeros(t, n, n)][nl.inputs, :]
    rhs2 = zeros(t, length(nl.inputs))

    # set the net currents to zero,
    eq3 = [zeros(t, n, n) I]
    rhs3 = zeros(t, n, 1)

    # but remove equations where we set the voltage and
    # thus don't know the net current
    i = setdiff(1:n, nl.inputs)
    eq3 = eq3[i,:]
    rhs3 = rhs3[i]

    # and finally, all injected currents (including those from
    # voltage inputs must balance to zero
    eq4 = [zeros(t, n)' ones(n)']
    rhs4 = zeros(t, 1)

    eq = [eq1; eq2; eq3; eq4]
    rhs = [rhs1; rhs2; rhs3; rhs4]

    (frequency, inputs, parameters) -> begin
        # fill in the admittances
        let ω = 2π * frequency * -1im
            for (coords, link) in nl.links
                r, l, c = resolve.([link.r, link.l, link.c], Ref(parameters))
                z = r + ω * c + l / ω
                # fill lower triangle
                eq[coords...] = z
                # fill upper triangle by reversing each (i,j) pair
                eq[reverse(coords)...] = z
            end

            # set the diagonal values to the negative sum of the other values
            i = diagind(eq)[1:n]
            eq[i] .= 0
            eq[i] = -sum(eq[1:n, 1:n], dims=1)
        end

        # and set the inputs
        m = length(nl.inputs)
        rhs[n+1:n+m] = inputs
        eq \ rhs
    end
end


"""
Given a `Netlist` and the parameters defining component values in that
circuit, `solve` returns a function that simulates a circuit. The
arguments of the returned function are frequency and the values of the
circuit parameters. The returned value is a vector of the output voltages.
"""
function solve(nl::Netlist, parameters::AbstractVector=[])
    s = raw(nl)
    (frequency, inputs) -> begin
        vi = s(frequency, inputs, parameters)
        return getindex.(Ref(vi), nl.outputs)
    end
end

"""
Returns a function that takes a frequency and computes the impedance
of a circuit between the two inputs. It is assumed that there are two
inputs, but ignores the outputs.
"""
function inputImpedance(nl::Netlist, parameters::AbstractVector=[])
    s = solve(nl, parameters)
    frequency -> begin
        vi = s(frequency, [1, 0])
        return 1/vi[nl.n+1]
    end
end

"""
Returns a function that computes the voltage transfer function of a
ciruit with an input, a reference voltage and an output at a
particular frequency.
"""
function transfer(nl::Netlist, parameters::AbstractVector=[])
    s = solve(nl, parameters)
    frequency -> begin
        return only(s(frequency, [1,0]))
    end
end

"""
Parses a numerical value with SI multiplier (and variants like u instead of μ)
"""
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



# A collection of simple functions useful in electronic design

"Inductance in nH of an air-core coil with n turns d mm in diameter and l mm long"
coil(n, d, l) = (l * n)^2 / (0.44 * d + l)

"""
Returns the inductance in nH for a specified core and number of turns. Throws
`ErrorException` if the core is unknown.
"""
function toroid(core::AbstractString, turns)
    if !(core in keys(toroid_parameters))
        error("Unknown toroid: $core")
    end
    toroid_parameters[core] * turns^2
end

# Toroid parameters after https://toroids.info (see https://kitsandparts.com/
# for compatible parts)
toroid_parameters = Dict(
    "T25-2"=> 3.4, "T25-6"=> 2.7, "T30-2"=> 4.3, "T30-6"=> 3.6,
    "T30-10"=> 2.5, "T37-0"=> 0.49, "T37-1"=> 8, "T37-2"=> 4,
    "T37-6"=> 3, "T37-7"=> 3.2, "T37-10"=> 2.5, "T37-17"=> 1.5,
    "T44-2"=> 5.2, "T44-6"=> 4.2, "T50-1"=> 10, "T50-2"=> 4.9,
    "T50-3"=> 17.5, "T50-6"=> 4, "T50-7"=> 4.3, "T50-10"=> 3.1,
    "T50-17"=> 1.8, "T68-1"=> 11.5, "T68-2"=> 5.7, "T68-6"=> 4.7,
    "T68-7"=> 5.2, "T68-10"=> 3.2)


end

# MicroSpice

MicroSpice is a very simple circuit simulator that simulates an electronic
circuit composed a variety of passive components. This circuit can be
constructed using a subset of standard Spice netlist syntax.

As an example, here is a simple RLC circuit with `in` as input and `out` as
output

```
L1 in  out 100nH
R1 out gnd  50Ω
C1 out gnd 100nF
```

Currently the only supported components are resistors (R), inductors (L)
and capacitors (C). Component values can be specified using standard SI
prefixes such as `k` (kilo = 1e3), `M` (Mega=1e6), `G`(giga=1e9), `m`
(milli=1e-3), `μ` (micro=1e-6), `n` (nano=1e-9), or `p` (pico=1e-12). In
addition, the non-standard forms `u` (micro=1e-6) and `meg` (Mega=1e6) are
used for Spice compatibility.

Units including H (Henry), F (Farad) and Ω (ohm) are allowed, but are
optional, but the unit must match the type of component.

# Running a Simulation

A `Netlist` describes an abstract circuit. To simulate that circuit,
the circuit is turned into a matrix equation consisting of complex
admittances that relate the voltages and currents in the circuit. In this
conversion, the voltages for any subset of nodes can be forced and currents
can be injected into any nodes that are not forced to a particular voltage.
The result of the voltages and all of the node-to-node currents are then
computed.

In normal practice, the internals of this conversion are hidden from view.
All you need to do is provide a `Netlist` and a list of inputs and outputs
and you get a simulation function that takes a frequency, a vector of
voltages and returns a vector of output voltages. You can call that
simulation function repeatedly for different frequencies or different
inputs.

As an example, here is a simple RLC circuit in Spice netlist format

```
L1 in  out 100nH
R1 out gnd  50Ω
C1 out gnd 100nF
```

If we take node `in` as the input, node `out` as the output and `gnd` as
the ground reference, we can simulate the voltage transfer function for
this circuit with this code:

```jldoctest; filter = r"(\d*)\.(\d{9})\d+" => s"\1.\2***"
nl = MicroSpice.Netlist("L1 in  out 100nH\nR1 out gnd 50Ω\nC1 out gnd 100nF\n")
s = MicroSpice.solve(nl, [:in => 1, :gnd => 0], "out")
decibel(x) = 20 * log10(abs(x))
(decibel ∘ s).([1.4e6, 1.5e6, 1.62e6, 1.8e6])
# output
4-element Vector{Float64}:
 12.883077832402897
 18.914298710446417
 27.65586909592151
 11.05634871505561
```

# Simulation Method

Once a `Netlist` is instantiated as a circuit consisting of ``n`` nodes
connected by complex admittances (the inverse of impedances), we can build
a system of linear equations in terms of the node voltages and the currents
injected into each node from the outside. Roughly half of this system
expresses the constraints of the circuit itself while the other half
expresses the voltages or currents imposed on the circuit.

The variables in the system are the node voltages and the injected
currents. For each internal node without any imposed voltage, we have a
current balance condition between internal node-to-node currents and the
(possibly zero) injected current. For nodes with an imposed voltage, we
don't know the injected current, but we do have a constraint on the node
voltages.

To make the system solvable, the sum of all injected currents must be zero
and the voltage of at least one node in each connected sub-circuit must be
set to a known value.

These linear equations are defined in several batches in terms of the ``n``
node voltages, the ``0 < k \le n`` forced voltages and the ``n-k``
known injected currents (any number of which can be zero).

## Current balance for nodes

The first batch of equations use Ohm's law to relate the voltages at each
node to the currents between these nodes. The sum of these node-to-node
currents is equated to the injected current.

The current ``i_{ab} `` from node ``a`` to node ``b`` is
```math
i_{ab} = \frac {v_a - v_b} {x_{ab}}
```
For all nodes that are not directly connected, ``i_{ab} = 0``. 

In addition, the current ``i_a`` injected into node ``a`` and the total leaving that node to other nodes will be equal. Thus,

```math
I_a = \sum_b i_{ab}
```

This can be rearranged and rewritten using admittances ``\sigma_{ab}``
instead of impedances ``x_{ab}`` into matrix form

```math
\left[
\begin{array}{ccccc}
 -\sum_{b\ne 1} \sigma_{1b} & \sigma_{12}&\ldots& \sigma_{1n}  \\ 
\sigma_{21} &-\sum_{b\ne 2} \sigma_{2b}&  \ldots& \sigma_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
\sigma_{n1} & \sigma_{n2}  & \ldots  & -\sum_{b\ne n} \sigma_{nb}
\end{array}
\;\middle|\;
-I
\right]
\begin{bmatrix}
v_1 \\
\vdots \\
v_n \\
i_1 \\
\vdots \\
i_n
\end{bmatrix} = 0 
```

## Nodes with forced voltages

Some possibly empty set of nodes will have voltages forced to a particular
value. This is expressed as

```math
v_a = V_a
```

In matrix form, each equation like this involves a row of coefficients with
a single non-zero value of ``1`` in the ``a``-th position.

At least node should be forced to a known voltage per connected
sub-circuit. In the extreme case, all of the nodes can be forced to known
voltages, although that isn't very interesting.

Note that forcing a node to a particular voltage involves injecting current
into the circuit at that point. By convention in MicroSpice, the currents
injected into nodes that are forced to a known voltage are unknown.
Conversely, the voltages for nodes with known injected currents are also
unknown.

## Current balance

The injected current for all nodes other than those with imposed voltages are set to zero (no injection) or to the known value. The equations for these current constraints are simple:

```math
i_a = I_a
```

As with the forced voltages, this is expressed in matrix form as a row per
constraint with a single non-zero coefficient equal to ``1``.

Note that, for simplicity, we assume that known current injection only
happens for nodes that are not forced to a particular voltage. This is not
supported in MicroSpice to avoid having to have an interface that allows
for current injection of unknown magnitude for a node with unknown voltage.

There is a single additional current balance that expresses the fact that all injected currents must balance. This is expressed with

```math
\sum_a i_a = 0
```
## Simple example

Taking the simple RLC circuit given above, the equations of state for
the circuit at a frequency ``f`` would be
```math
\begin{bmatrix}
-\frac j {2 \pi f L} & \frac j {2 \pi f L}& 0 & -1 & 0 & 0\\
\frac j {2 \pi f L} & -\frac j {2 \pi f L} -\frac 1 R - 2\pi jf C &\frac 1 R + 2\pi jf C & 0 & -1 & 0\\
0 & \frac 1 R + 2 \pi jf C & -\frac 1 R - 2\pi jf C & 0 & 0 & -1\\ 
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\ 
0 & 0 & 0 & 0 & 1 & 0 \\ 
0 & 0 & 0 & 1 & 1 & 1 \\
\end{bmatrix}
\begin{bmatrix}
v_{in} \\ v_{out} \\ v_{gnd} \\
I_{in} \\ I_{out} \\ I_{gnd}
\end{bmatrix}
=
\begin{bmatrix}
0 \\ 0 \\ 0 \\ 1 \\ 0 \\ 0 \\ 0
\end{bmatrix}
```
Where ``R = 50\Omega``, ``C = 100nF``, and ``L = 100nH``.

## Parametrized circuits

The value for a component in a `Netlist` can also be a reference to a
named parameter such as `$R3`. The names for such parameters are
passed in when a netlist is created and the values are passed in when
the circuit is materialized. At that point, you can solve the circuit
for different frequencies.

Parametrizing a circuit is handy, for instance, if you want to investigate
the effect of component variation or if you want to use an optimizer
to find optimal component values for complex objectives. Importantly,
the parameters are injected directly from your values into the circuit
which should allow auto-differentiation during optimization.

As an example, let's take the example above, but analyze it for the effect
of ±5% variation of the nominal component values of 100nH and 100nF.
```julia
nl = MicroSpice.Netlist(raw"""
L1 in  out $L
R1 out gnd 50Ω
C1 out gnd $C
""", ["L", "C"])
fx = 0.5e6:0.005e6:2.5e6
nominal = [
	decibel(only(MicroSpice.solve(nl, [:in, :gnd], ["out"], [100e-9, 100e-9])(f, [1, 0])))
	for f in fx
	]
r = [ 
       decibel(only(MicroSpice.solve(nl, [:in, :gnd], ["out"], Vector(params))(f, [1, 0])))
       for f in fx,
       params in  eachrow((0.95 .+ 0.1 * rand(50,2)) .* [100e-9, 100e-9]')
       ]
```
At this point, `r` has an array with 50 columns, one for each
combination of parameter values. We can plot the nominal response
surrounded by the range of possibilities like this:
```julia
x = vcat(fx./1e6, reverse(fx./1e6))
y = vcat(vec(minimum(r, dims=2)),reverse(vec(maximum(r, dims=2))))
plot(Shape(x,y), fill="lightblue", width=0, ylab="Voltage gain (db)",xlab="Frequency (MHz)")
plot!(fx./1e6, nominal, legend=false)
```
Note that most of the variation in response is due to shifts in the
resonant frequency, not changes in the peak response.

The resulting graph looks like this:

<img width="592" alt="image" src="https://github.com/user-attachments/assets/23c048bc-422c-4e18-824b-af75cc09438c">


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tdunning.github.io/MicroSpice.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tdunning.github.io/MicroSpice.jl/dev/)
[![Build Status](https://github.com/tdunning/MicroSpice.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tdunning/MicroSpice.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tdunning/MicroSpice.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tdunning/MicroSpice.jl)

using MicroSpice
using Plots
using Statistics

evalCircuit(nl, frequencies, penalties) = parameters -> quality(nl, parameters, frequencies, penalties)

"""
Finds a worst case value of performance over variation of component values
"""
function quality(sim, parameters, frequencies, penalties)
    let n = 100
        maximum(i -> begin
                    p = parameters .* (1 .+ 0.05 *randn(size(parameters)))
                    q0(sim, p, frequencies, penalties)
                end, 1:n)
    end
end

"""
Internal function that returns the fitness of a solution prioritizing
the net attenuation at frequencies of interest adjusted by penalties.
There is a secondary goal the encourages low attenuation at the primary
frequency.
"""
function q0(sim, parameters, frequencies, penalties)
    response = [20 * log10(abs(only(sim(f, [2,0], parameters)))) for f in frequencies]
    net =  (response[2:end] .- response[1]) .+ penalties
    return maximum(net) - response[1]
end

"""
State of an evolutionary search species
"""
struct ClimbState
    ix::Vector         # the current selection
    old::Vector        # the previous (for recorded step)
    value::Float64     # performance at current point
    rate::Float64      # omni mutation rate
end

hash(x::ClimbState) = hash(x.ix)

"""
Mutates a current state into a new state given

c - the old state
q - the fitness function (lower is better)
parameterSets - a list of lists of parameter values to translate indexes in c into values for q
μ - how much history should affect the new mutation rate (near one for lots of history, near zero for none)

"""
function mutate(c::ClimbState, q, parameterSets, μ)::ClimbState
    ix = copy(c.ix)
    rate = c.rate
    omni = 0
    m = length(parameterSets)
    for i in 1:m
        # this is like exponential distribution but with a long tail
        v = -log(1-rand()) + log(1-rand())^2
        omni += v
        # find new state with the omni mutation plus a move in the previous direction(ish)
        ix[i] += Int(round(rand(-1:2:1) * rate * v + (1 + randn()) * (ix[i] - c.old[i])))
        # avoid running off the ends
        if ix[i] < 1
            ix[i] = 1
        elseif ix[i] > length(parameterSets[i])
            ix[i] = length(parameterSets[i])
        end
    end
    # the new rate is a combination of historical rates and the latest jump
    # the recorded step also feeds into the new omni mutation
    rate = μ * rate + (1-μ) * 0.8 * (rate * omni / m + mean(ix .- c.ix) / 10)
    return ClimbState(ix, c.ix, q(getindex.(parameterSets, ix)), rate)
end

"""
This is an meta-evolutionary hill-climber. Each climber has a current
state and an evolution rate. To mutate a climber, we combine omni-directional
mutation scaled by the rate component of a mutation state plus some 
mutation along the direction of the previous step. The size of the 
omni-directional mutation is then used to mutate the scale of the next
omni-directional mutation. In each generation
the fraction of survivors is λ.

The parameters are:

q - a function that returns a quality function 
ix - an array of indexes into the parameterSets that define the starting point
parameterSets - an array of arrays of component values. This translates ix into what q wants 

Optional keyword parameters include:

pop - the population size
λ - the survival fraction
μ - momentum for omni-directional mutation rate, 1 means stable, 0 means latest value rules
rate - initial value for the omni-directional mutation
gens - the number of generations 
info - how often to report top few values (small means often)
"""
function metaClimb(q, ix, parameterSets; pop=100, λ=0.2, μ=0.8, rate=1, gens=30, info=20)
    if info < 1
        info = 1
    end
    history = []

    # pre-allocate simulation functions to allow multi-threading
    sims = [q() for i in 1:Threads.nthreads()]

    # fill up the population by mutating the intial state
    population = Vector{ClimbState}()
    c0 = ClimbState(ix, ix, 0.0, rate)
    while length(population) < pop
        push!(population, mutate(c0::ClimbState, sims[1], parameterSets, μ))
    end
    
    survivors = Int(floor(λ * pop) + 1)
    j = 1
    for g in 1:gens
        if g%info == 0
            @info "step" g population[1] population[2]
        end
        # find the best examples
        sort!(population, by=c->c.value + c.rate/100)
        push!(history, population[1:10])
        # apart from the survivors, fill out with new mutants
        Threads.@threads for i in survivors+1:pop
            c = mutate(population[j], sims[Threads.threadid()], parameterSets, μ)
            population[i] = c
            j = j % survivors + 1
        end
    end
    # sort one last time to show the actual results
    sort!(population, by=c->c.value + c.rate/100)
    return history, population
end


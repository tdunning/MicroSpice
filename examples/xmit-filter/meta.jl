using LinearAlgebra

"""
Recorded step meta evolution as described in

https://arxiv.org/abs/0803.3838

The basic idea is that we have a cloud of candidate solutions.  At
each step, we take the best few of these, mutate them and rank
according to a quality function.

The trick with recorded step evolution is that each candidate only
only has a potential optimum for our quality function, it also has a
representation of a mutation rate. The algorithm learns about good
solutions at the same time as it learns about good mutation rates (and
directions). This improves convergence over evolutionary approaches
without meta-mutation.
"""

struct Candidate{T} 
    omni::Float64
    state::Vector{T}
    old::Vector{T}
end

struct MetaEp
    population::Vector{Candidate{Float64}}
    len::Integer
    retention::Float64
    fitness::Function
    offset::Float64
    
    function MetaEp(x0::Matrix, retention, fitness; len=size(x0,1), omni=1.0, offset=1)
        p = Vector{Candidate{Float64}}()
        for row in eachrow(x0)
            c = Candidate(omni, copy(row), copy(row))
            append!(p, Ref(c))
        end
        extend(p, len, omni)
        return new(p, len, retention, fitness, offset)
    end
end

function step(pop::MetaEp)
    r = [pop.fitness(z.state) for z in pop.population]
    k = sortperm(r)
    n = Int(pop.retention * length(pop.population))
    p = copy(pop.population[k[1:n]])
    extend(p, pop.len, pop.offset)
    pop.population .= p
end

function extend(p::Vector, len, offset)
    i = 1
    base = length(p)
    while length(p) < len
        i = i%base + 1
        append!(p, Ref(mutate(p[i], offset)))
    end
end

function mutate(c::Candidate, offset)::Candidate
    step = c.state - c.old
    next = c.state .+ c.omni .* randn(length(c.state)) .+ (randn() + offset) .* step
    return Candidate(-1.5 * (c.omni + 0.05*norm(next-c.state)) * log(1-rand()), next, c.state)
end

function bowl(center::Vector)::Function
    v -> let d = (v .- center)
        sum(d .* d)
    end
end

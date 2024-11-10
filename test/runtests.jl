using MicroSpice
using Test, LinearAlgebra


@testset "decoding values" begin
    @test MicroSpice.decode('R', "100") ≈ MicroSpice.Res(100)
    @test MicroSpice.decode('R', "1e5") ≈ MicroSpice.Res(100_000)
    @test MicroSpice.decode('R', "100k") ≈ MicroSpice.Res(100e3)
    @test MicroSpice.decode('R', "100μ") ≈ MicroSpice.Res(100e-6)
    @test MicroSpice.decode('R', "100 ohm") ≈ MicroSpice.Res(100)
    @test MicroSpice.decode('R', "100 G ohm") ≈ MicroSpice.Res(100e9)
    @test MicroSpice.decode('R', "100megaΩ") ≈ MicroSpice.Res(100e6)
    @test MicroSpice.decode('R', "100mega") ≈ MicroSpice.Res(100e6)
    @test MicroSpice.decode('C', "10nF") ≈ MicroSpice.Cap(10e-9)
    @test MicroSpice.decode('C', "10p") ≈ MicroSpice.Cap(10e-12)
    @test MicroSpice.decode('C', "10MF") ≈ MicroSpice.Cap(10e6)
    @test MicroSpice.decode('C', "10megF") ≈ MicroSpice.Cap(10e6)
    @test MicroSpice.decode('L', "100pH") ≈ MicroSpice.Ind(100e-12)
    @test MicroSpice.decode('L', "100n") ≈ MicroSpice.Ind(100e-9)

    @test_throws ErrorException MicroSpice.decode('R', "100F")
    @test_throws ErrorException MicroSpice.decode('R', "1H")
    @test_throws ErrorException MicroSpice.decode('C', "2e3Ω")
    @test_throws ErrorException MicroSpice.decode('C', "2e3 ohm")
    @test_throws ErrorException MicroSpice.decode('L', "5F")
    @test_throws ErrorException MicroSpice.decode('L', "2e3Ω")
    @test_throws ErrorException MicroSpice.decode('L', "2e3 ohm")
    @test_throws ErrorException MicroSpice.decode('X', "2e3")
end

@testset "internal expressions" begin
    @test MicroSpice.resolve(MicroSpice.InternalExpression('R', 3, []), []) ≈ 3
    @test MicroSpice.resolve(MicroSpice.InternalExpression('R', 3), []) ≈ 3
    @test MicroSpice.resolve(MicroSpice.InternalExpression('C', 3, []), []) ≈ 3
    @test MicroSpice.resolve(MicroSpice.InternalExpression('C', 3), []) ≈ 3
    @test MicroSpice.resolve(MicroSpice.InternalExpression('L', 3, []), []) ≈ 3
    @test MicroSpice.resolve(MicroSpice.InternalExpression('L', 3), []) ≈ 3
    @test MicroSpice.resolve(MicroSpice.Res(3).r, []) ≈ 1/3
    @test MicroSpice.resolve(MicroSpice.Cap(3).c, []) ≈ 3
    @test MicroSpice.resolve(MicroSpice.Ind(3).l, []) ≈ 1/3
    let r1 = MicroSpice.InternalExpression('R', :r1, ["a", "r1", "r2", "l", "c"]),
        r2 = MicroSpice.InternalExpression('R', "r2", [:l1, :r1, :r2, :l2, :c1, :c2]),
        l1 = MicroSpice.InternalExpression('L', "l1", [:l1, :r1, :r2, :l2, :c1, :c2]),
        l2 = MicroSpice.InternalExpression('L', "l2", [:l1, :r1, :r2, :l2, :c1, :c2]),
        c1 = MicroSpice.InternalExpression('C', "c1", [:l1, :r1, :r2, :l2, :c1, :c2]),
        c2 = MicroSpice.InternalExpression('C', "c2", [:l1, :r1, :r2, :l2, :c1, :c2])

        @test_throws BoundsError MicroSpice.resolve(r1, []) ≈ 1/34
        @test MicroSpice.resolve(r1, [rand(), 34, rand(4)...]) ≈ 1/34
        @test MicroSpice.resolve(r1 + r2, [rand(), 2, 3, rand(3)...]) ≈ 1/2 + 1/3
        @test MicroSpice.resolve(c1 + c2, [rand(4)..., 6, 8]) ≈ 14.0
        @test MicroSpice.resolve(l1 + l2, [10, rand(2)..., 5, rand(2)...]) ≈ 1/10 + 1/5
    end
end
        

@testset "construction of a vanilla netlist" begin
    let nl = MicroSpice.Netlist("L1 in  out 100n\nR1 out gnd 50\nC1 out gnd 100n\n")
        @test nl.n == 3
        @test getindex.(Ref(nl.names), ["in", "out", "gnd"]) == [1,2,3]

        @test nl.links[1,3] ≈ MicroSpice.Link(0, 0, 0)
        @test nl.links[3,1] ≈ MicroSpice.Link(0, 0, 0)
        for i in 1:3
            @test nl.links[i, i] ≈ MicroSpice.Link(0, 0, 0)
        end

        @test nl.links[2,3] ≈ MicroSpice.Link(0.02, 1.0e-7, 0)
        s = MicroSpice.solve(nl, [:in, :gnd], [:out])
        y = [20*log10(abs(only(s(f, [1, 0])))) for f in 1.2e6:0.1e6:1.8e6]
        @test y ≈ [7.294866, 9.545497, 12.883077, 18.914298, 32.859821, 16.921493, 11.056348] atol=1e-5
        x = MicroSpice.inputImpedance(nl, :in, :gnd)
        @test abs.(x.([1.2e6, 1.5e6])) ≈ [0.5724569, 0.1202032] atol=1e-6
        x = MicroSpice.transfer(nl, :in, :out, :gnd)
        @test abs.(x.([1.2e6, 1.5e6])) ≈ [2.31602539, 8.82500450] atol=1e-6
    end
end


@testset "construction of a netlist with parameters" begin
    let nl = MicroSpice.Netlist("L1 in  out \$L\nR1 out gnd 50\nC1 out gnd 100nF\n", ["L"])
        @test_throws BoundsError MicroSpice.ComplexCircuit(nl, 1.5e6, [])
        c = MicroSpice.ComplexCircuit(nl, 1.5e6, [100e-9])

        # without parameters, the sim is doomed but the problem only shows at bind-time
        s = MicroSpice.solve(nl, [:in, :gnd], [:out])
        @test_throws BoundsError s(1e6, [1, 0])

        s = MicroSpice.solve(nl, [:in, :gnd], [:out], [100e-9])
        y = [20*log10(abs(only(s(f, [1, 0])))) for f in 1.2e6:0.1e6:1.8e6]
        @test y ≈ [7.294866, 9.545497, 12.883077, 18.914298, 32.859821, 16.921493, 11.056348] atol=1e-5

        s = MicroSpice.solve(nl, [:in, :gnd], [:out], [120e-9])
        y = [20*log10(abs(only(s(f, [1, 0])))) for f in 1.2e6:0.1e6:1.8e6]
        @test y ≈ [9.942511, 13.964684, 22.554476, 23.136562, 13.386021, 8.635949, 5.423041] atol=1e-5
    end
end

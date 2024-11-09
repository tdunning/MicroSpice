using MicroSpice
using Test, LinearAlgebra

@testset "construction of a netlist" begin
    let nl = MicroSpice.Netlist("L1 in  out 100n\nR1 out gnd 50\nC1 out gnd 100n\n")
        @test nl.n == 3
        @test getindex.(Ref(nl.names), ["in", "out", "gnd"]) == [1,2,3]
        @test sum(diag(nl.links)) + nl.links[3,1] + nl.links[1,3] == MicroSpice.Link(0,0,0)
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



@testset "decoding values" begin
    @test MicroSpice.decode('R', "100") ≈ 100
    @test MicroSpice.decode('R', "1e5") ≈ 100_000
    @test MicroSpice.decode('R', "100k") ≈ 100e3
    @test MicroSpice.decode('R', "100μ") ≈ 100e-6
    @test MicroSpice.decode('R', "100 ohm") ≈ 100
    @test MicroSpice.decode('R', "100 G ohm") ≈ 100e9
    @test MicroSpice.decode('R', "100megaΩ") ≈ 100e6
    @test MicroSpice.decode('R', "100mega") ≈ 100e6
    @test MicroSpice.decode('C', "10nF") ≈ 10e-9
    @test MicroSpice.decode('C', "10p") ≈ 10e-12
    @test MicroSpice.decode('C', "10MF") ≈ 10e6
    @test MicroSpice.decode('C', "10megF") ≈ 10e6
    @test MicroSpice.decode('L', "100pH") ≈ 100e-12
    @test MicroSpice.decode('L', "100n") ≈ 100e-9

    @test_throws ErrorException MicroSpice.decode('R', "100F")
    @test_throws ErrorException MicroSpice.decode('R', "1H")
    @test_throws ErrorException MicroSpice.decode('C', "2e3Ω")
    @test_throws ErrorException MicroSpice.decode('C', "2e3 ohm")
    @test_throws ErrorException MicroSpice.decode('L', "5F")
    @test_throws ErrorException MicroSpice.decode('L', "2e3Ω")
    @test_throws ErrorException MicroSpice.decode('L', "2e3 ohm")
end

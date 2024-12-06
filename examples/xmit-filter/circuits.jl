using MicroSpice

"""
Parallel LC resonant narrow-pass filter
"""
nl0 = MicroSpice.Netlist(raw"""
R1 in out 50
L1 out gnd $L1
C1 out gnd $C1
R2 out gnd 50
""", [:L1, :C1], [:in, :gnd], [:out])
c0 = [40e-9, 800e-12]

"""
Low-pass π-filter with RLC on output shunt
leg
"""
nl1 = MicroSpice.Netlist(raw"""
R1 in N001 50
L1 N001 out $L1
C1 N001 gnd $C1
C2 out N002 $C2
R2 out gnd 50
L2 N002 gnd $L2
R3 N002 gnd $R
""", [:L1, :C1, :C2, :L2, :R], [:in, :gnd], [:out])

"""
Low-pass π-filter
"""
nl2 = MicroSpice.Netlist(raw"""
R1 in N001 50
L1 N001 out $L1
C1 N001 gnd $C1
C2 out gnd $C2
R2 out gnd 50
""", [:L1, :C1, :C2], [:in, :gnd], [:out])
c2 = [300e-9, 150e-12, 150e-12]

"""
LC series shunt
"""
nl3 = MicroSpice.Netlist(raw"""
R1 in out 50
C1 out N001 225p
R2 out gnd 50
L2 N001 gnd 20n
R3 N001 gnd 30
""", [], [:in, :gnd], [:out])

"""
Chebyshev 5th order with 5db passband ripple
formed with double-π
"""
nl4 = MicroSpice.Netlist(raw"""
R1 in N001 50
L1 N001 N002 $L1
C2 N001 gnd $C1
C1 N002 gnd $C2
L2 N002 out $L2
C3 out gnd $C3
R2 out gnd 50
""", [:L1, :L2, :C1, :C2, :C3], [:in, :gnd], [:out])
c4 = [162.1753e-9, 162.1753e-9, 550.6910e-12, 700.1618e-12, 550.6910e-12]

"""
Chebyshev 3rd order with 5dB passband ripple formed with 
π-filter
"""
nl5 = MicroSpice.Netlist(raw"""
R1 in N001 50
C1 N001 gnd $C1
L1 N001 out $L1
C2 out gnd $C2
R2 out gnd 50
""", [:L1, :C1, :C2], [:in, :gnd], [:out])
c5 = [151.9541e-9, 529.6833e-12, 529.6833e-12]

"""
Inverse Chebyshev 3rd order formed as π with series element
formed by parallel LC element
"""
nl6 = MicroSpice.Netlist(raw"""
R1 in N001 100
L1 N001 out $L1
C1 N001 gnd $C1
C2 N001 out $C2
C3 out gnd $C3
R2 out gnd 200
""", [:L1, :C1, :C2, :C3], [:in, :gnd], [:out])
c6 = [500e-9, 100e-12, 10e-12, 100e-12]


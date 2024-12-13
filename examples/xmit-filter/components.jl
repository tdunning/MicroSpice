using MicroSpice

# inventory of my parts box
# toroid inductance formula from https://toroids.info/T25-6.php and https://toroids.info/T37-6.php
inductors(includeHandMade=false) = sort(collect(Set([
    # discrete SMD components
    1.8e-9, 2.2e-9, 2.7e-9, 3.3e-9, 4.7e-9, 5.6e-9, 6.8e-9,
    8.2e-9, 10e-9, 15e-9, 18e-9, 27e-9,
    5.5e-9, 6.0e-9, 8.9e-9, 12.3e-9, 15.7e-9, 19.4e-9,
    6.9e-9, 10.2e-9, 11.2e-9, 13.7e-9, 17.0e-9, 22e-9,
    8.1e-9, 12.1e-9, 14.7e-9, 16.6e-9, 21.5e-9, 23e-9,
    25e-9, 27.3e-9,
    330e-9, 560e-9, 820e-9, 1.2e-6, 2.2e-6, 4.7e-6, 6.8e-6,
    if includeHandMade
        [
            # T37-6 toroid with 2 to 12 turns
            1e-9 * [MicroSpice.toroid("T37-6", t) for t in 2:12]...,
            # T25-6 toroid with 2 to 9 turns
            1e-9 * [MicroSpice.toroid("T25-6", t) for t in 2:12]...,
            # wound inductors on chopstick assuming last turn only goes halfway
            1e-9 * [MicroSpice.coil(n+0.5, 5, 0.6*n) for n in 3:12]...
                ]
    else
        []
    end...
        ])))

caps() = sort(collect(Set(
[
    # SMD capacitors
    15e-12, 22e-12, 27e-12, 33e-12, 47e-12, 68e-12,
    82e-12, 100e-12, 120e-12, 150e-12, 180e-12,
    1e-12, 1.5e-12, 1.8e-12, 2.2e-12, 2.7e-12, 3.3e-12,
    3.9e-12, 4.7e-12, 5.6e-12, 6.8e-12, 8.2e-12, 10e-12
])))

resistors() = sort(collect(Set(
    [
        # SMD resistors
        1, 1.2, 1.5, 2, 2.7, 3.3, 4.3, 5.1, 6.8, 8.2, 10, 12, 15, 20, 27,
        33, 43, 51, 68, 82, 100, 120, 150, 200, 270, 330, 40, 510, 680,
        820, 1e3, 1.2e3, 1.5e3, 2e3, 2.7e3, 3.3e3, 4.3e3, 5.1e3, 6.8e3,
        8.2e3, 10e3
    ])))

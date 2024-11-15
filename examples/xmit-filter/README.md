# Transmit filtering for direct frequency synthesis

Chips like the Si5351 all very simple systems to directly generate
frequencies up into the VHF band or even into UHF if the third
harmonic is used. Unfortunately, the generated signal is (roughly) a
square wave which means that before you connect it to an antenna, you
need to filter out undesired harmonics.

There are excellent design tools for that help you calculate a nice
filter, but what happens when you need to specialize to the particular
set of parts in your parts box? What happens when those components
vary with temperature or time? Is there a clever design hack that will
suppress harmonics more than a traditional filter might?

All of these questions and more can be answered with MicroSpice.

As a specific example, let's design a filter to eliminate harmonics
from a WSPR transmission in the 10m band (28.1261 MHz ± 100Hz). As a
wrinkle to this problem, we want to account for parts variations and
non-ideal input waveforms.

# What does success mean?

If you look at the spectrum of an idealized square wave, the second
harmonic is completely suppressed and the third and fifth harmonics
are roughly 9.5 and 14 dB below the fundamental. If the duty cycle of
the square wave is not exactly 50%, however, the second harmonic may
become significant. What this means is that we want to suppress
different harmonics by different amounts because they will be, to
differing extents already suppressed. For instance, if we look at the
[spectrum of a 55% duty cycle waveform](https://gist.github.com/tdunning/e551ae973422609f031c1bdca39ff5b4)
with 1ns rise and fall times, we see that the 3rd and 5th harmonics
are at roughly -10dBc and -15dBc as was the case with the ideal square
wave, but the 2nd and 4th harmonics are now non-negligible at just
above -20dBc (with the 4th a few dB _higher_ than all but the 3rd).

<img width="400" alt="A plot of the spectrum of a 55% duty cycle square wave with 1ns rise and fall times" src="https://gist.github.com/user-attachments/assets/41428c2b-e634-4fb3-b912-e1b0a97428cd">

This is interesting because it raises the possibility that we can
design a filter that uses this property to advantage. That filter
might produce better results than we might expect or be simpler for
the same performance.

To do these fancy optimizations requires more than a canned filter
design. Software for filter design rarely allows for such nuanced
objectives, but it can be used to get in the ball park.

```spice
R1 in N001 50
L1 N001 out $L1
C1 N001 gnd $C1
C2 out N002 $C2
R2 out gnd 50
L2 N002 gnd $L2
R3 N002 gnd $R
```

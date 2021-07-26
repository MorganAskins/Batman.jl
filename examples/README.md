Batman.jl Examples
==================

### 0_Batsignal
Using the wolfram batsignal function, this example produces a random sample
of the function with the addition of noise.

### 1_Poisson
Simple example of a counting experiment where each of the parameters includes
a fixed uncertainty on the rate.

### 2_SpectralFit
Similar to the counting example, except that each parameter has an associated
shape to further distinguish it from the others.

### 3_UpRoot
Using the [UpRoot] library, [ROOT] files are read in and compared to similar
distributions created directly in Julia. An extended likelihood fit to the
shapes of the spectra is preformed.

### 4_SpectralMonofit


[UpRoot]: https://github.com/JuliaHEP/UpROOT.jl
[ROOT]: https://root.cern.sh

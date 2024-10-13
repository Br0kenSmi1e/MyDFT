using MyDFT

x = SpaceSample1D(-10.0,10.0,100)
sys = InteractingHamiltonian1D(x,ξ->ξ^2)

n = self_consistent_loop(sys,17,periodic=true)

using CairoMakie
lines(coordinate_sample(x),n)
using MyDFT
using MyDFT: NonInteractingHamiltonian, InteractingHamiltonian

T = 100
x = collect(range(-π,π,length=T))

h = NonInteractingHamiltonian(x->x^2, -π, 100, 2π/100)
num_electron = 16
n_non_inter = self_consistent_loop(h, num_electron)

hi = InteractingHamiltonian(x->x^2, -π, 100, 2π/100)
n_inter = self_consistent_loop(hi, num_electron)

using CairoMakie
f = Figure()
ax = Axis(f[1,1])
lines!(ax,x,n_non_inter,color = :red)
lines!(ax,x,n_inter,color = :black)
f
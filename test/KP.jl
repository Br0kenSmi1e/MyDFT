using MyDFT
using CairoMakie

x = SpaceSample1D(0.0, 1.0, 100)
sys = one_dim_kronig_penney(0.5, 1e6, x)

n, E, ψ = kohn_sham_solver(sys, 1, zero(coordinate_sample(x)))

f = Figure()
ax = Axis(f[1,1])
for i in 1:7
    lines!(ax, coordinate_sample(x), ψ[:,i])
end
f
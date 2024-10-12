using MyDFT
using LinearAlgebra
using CairoMakie

function harmonic_potential(x::T) where T<:Real
    return x*x/2
end

T = 100
x = collect(range(-2π,2π,length=T))
spec = eigen(Matrix(one_dim_kinetic_energy_open(T,4π/T))+Matrix(one_dim_external_potential(harmonic_potential,x)))

f = Figure()
ax = Axis(f[1,1])
for k in 1:5
    lines!(ax,x,spec.vectors[:,k])
end
f
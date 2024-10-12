using MyDFT
using CairoMakie
using LinearAlgebra
using KrylovKit

# one_dim_kinetic_energy(10,1)
T = 40
x = range(0,2π,length=T)
# plane_wave = sin.(x)
# kinetic_energy = one_dim_kinetic_energy(T,2π/T)*plane_wave
spec = eigen(Matrix(one_dim_kinetic_energy_periodic(T,2π/T)))


f = Figure()
ax = Axis(f[1,1])
# lines!(ax,x,plane_wave)
# lines!(ax,x,kinetic_energy)
for k in 1:5
    lines!(ax,x,spec.vectors[:,k])
end
f
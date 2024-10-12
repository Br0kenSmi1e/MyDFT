# this file deals with one dimensional electronic systems
# the functions are sampled on a grid {x₀,x₀+h,x₀+2h,...,x₀+Gh} (a G+1-dimesional vector)
# and the function f is represented by a vector {f(x₀),f(x₀+h),...,f(x₀+(G-1)h)}
# Also, we assume periodic boundary condition, f(x₀+Gh) = f(x₀)

abstract type Hamiltonian{T} end

function one_dim_kinetic_energy_periodic(G::Int, h::T; open::Bool=false) where T<:Real
    # f''(x₀+kh) = [f(x₀+(k-1)h) + f(x₀+(k+1)h) - 2f(x₀+kh)] / h² + O(h²), 1≤k≤G-1
    operator = Matrix(Tridiagonal(-0.5 .* ones(G-1)./(h^2), ones(G)./(h^2), -0.5 .* ones(G-1)./(h^2)))
    if !open
        operator[1,G] = -0.5/(h^2)
        operator[G,1] = -0.5/(h^2)
    end
    return operator
end

function one_dim_external_potential(V, x::AbstractVector{T}) where T<:Real
    # V is the external potential that this electronic system lives in
    return Diagonal([V(ξ) for ξ in x])
end

function one_dim_coulomb(n::Vector{T}, x::AbstractVector{T}, h::T) where T<:Real
    G = length(x)
    integrand = [[n[k]*h/sqrt((x[l]-x[k])^2+eps(T)) for k in 1:G] for l in 1:G]
    # @show sum(integrand,dims=1)
    return Diagonal(sum(integrand,dims=1)[1])
end

function one_dim_lda(n::AbstractVector{T}) where T<:Real
    # v(x) = -(3n(x)/π)^(1/3)
    return Diagonal(-(3 .* n ./ π).^(1/3))
end

function kohn_sham_solver(h::Hamiltonian{T}, density, num_electron::Int) where T<:Complex
    m = mat(h, density)
    @assert num_electron <= 2 * size(m, 1)
    @assert ishermitian(m)

    ψ = eigvecs(m)
    nl = zero(density)
    for k in 1:num_electron ÷ 2
        nl .+= 2 .* abs2.(view(ψ, :, k)) ./ step(h)
    end
    if num_electron % 2 == 1
        nl .+= abs2.(view(ψ, :, num_electron ÷ 2 + 1)) ./ step(h)
    end
    return nl
end

struct NonInteractingHamiltonian{FT, RT} <: Hamiltonian{Complex{RT}}
    potential::FT
    x0::RT
    num_coordinates::Int
    step::RT
end
struct InteractingHamiltonian{FT, RT} <: Hamiltonian{Complex{RT}}
    potential::FT
    x0::RT
    num_coordinates::Int
    step::RT
end
function coordinates(h::Hamiltonian)
    return LinRange(x0(h), x0(h) + num_coordinates(h) * step(h), num_coordinates(h))
end
for HT in [:InteractingHamiltonian, :NonInteractingHamiltonian]
    @eval num_coordinates(x::$HT) = x.num_coordinates
    @eval step(x::$HT) = x.step
    @eval x0(x::$HT) = x.x0
    @eval potential(x::$HT) = x.potential
end

function mat(h::NonInteractingHamiltonian{F, T}, density) where {F, T}
    return one_dim_kinetic_energy_periodic(num_coordinates(h), step(h)) +
        one_dim_external_potential(potential(h), coordinates(h))
end
function mat(h::InteractingHamiltonian{F, T}, density) where {F, T}
    x = coordinates(h)
    @show density
    return one_dim_kinetic_energy_periodic(h.num_coordinates, h.step) +
        one_dim_external_potential(h.potential, x)
        one_dim_coulomb(density, x, h.step) +
        one_dim_lda(density)
end


function self_consistent_loop(
        hamiltonian::Hamiltonian{T}, 
        num_electron::Int,
        tol = 1e-3, max_iter::Int = 10^4,
        n_ini = ones(real(T), num_coordinates(hamiltonian))
    ) where T<:Complex
    n_old = n_ini
    local n_new
    for _ in 1:max_iter
        n_new = kohn_sham_solver(hamiltonian, n_old, num_electron)
        sum(abs2, n_new - n_old) < tol && break
        n_old = n_new
    end

    return n_new
end
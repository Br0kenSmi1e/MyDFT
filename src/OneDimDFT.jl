struct SpaceSample1D{T}
    left_boundary::T
    right_boundary::T
    num_samples::Int
end

function step(X::SpaceSample1D)
    return (X.right_boundary - X.left_boundary) / (X.num_samples-1)
end

function coordinate_sample(X::SpaceSample1D)
    return collect(range(X.left_boundary, X.right_boundary, X.num_samples))
end

function kinetic_energy_1d(X::SpaceSample1D{T}; periodic::Bool = true) where T<:Real
    x_step = step(X)
    operator = Tridiagonal(ones(T, X.num_samples - 1) .* (-0.5/x_step^2),
        ones(T, X.num_samples) .* (1/x_step^2),
        ones(T, X.num_samples - 1) .* (-0.5/x_step^2))
    if periodic
        operator = Matrix(operator)
        operator[1, X.num_samples] = -0.5/x_step^2
        operator[X.num_samples, 1] = -0.5/x_step^2
    end
    return operator
end

function external_potential_1d(X::SpaceSample1D{T}, V_ext) where T<:Real
    return Diagonal([V_ext(x) for x in coordinate_sample(X)])
end

function coulomb_potential_1d(X::SpaceSample1D{T}, density::AbstractVector{T}) where T<:Real

    @assert length(density) == X.num_samples

    x = coordinate_sample(X)
    x_step = step(X)
    integrand = [[density[k]*x_step/sqrt((x[l]-x[k])^2+0.1) for l in 1:X.num_samples] for k in 1:X.num_samples]
    return Diagonal(sum(integrand))
end

function lda_potential_1d(density::AbstractVector{T}) where T<:Real
    return Diagonal((3/π)^(1/3) .* density .^(1/3))
end

abstract type Hamiltonian1D{T} end

struct NonInteractingHamiltonian1D{RT, FT} <: Hamiltonian1D{RT}
    space::SpaceSample1D{RT}
    external_potential::FT
end

struct InteractingHamiltonian1D{RT, FT} <: Hamiltonian1D{RT}
    space::SpaceSample1D{RT}
    external_potential::FT
end

function mat(system::NonInteractingHamiltonian1D{RT, FT},
        density::AbstractVector{RT}; periodic::Bool=true
        ) where {RT<:Real, FT}
    return kinetic_energy_1d(system.space,periodic=periodic) +
        external_potential_1d(system.space, system.external_potential)
end

function mat(system::InteractingHamiltonian1D{RT, FT},
        density::AbstractVector{RT}; periodic::Bool=true
        ) where{RT<:Real, FT}
    @assert length(density) == system.space.num_samples

    return kinetic_energy_1d(system.space,periodic=periodic) +
        external_potential_1d(system.space, system.external_potential) +
        coulomb_potential_1d(system.space, density) +
        lda_potential_1d(density)
end

function kohn_sham_solver(H::Hamiltonian1D{T},
        num_electron::Int,
        density::AbstractVector{T};
        periodic::Bool=true
        ) where T<:Real
    h = mat(H, density, periodic=periodic)
    @assert length(density) == H.space.num_samples
    @assert ishermitian(h)

    x_step = step(H.space)
    E,ψ = eigen(h)
    density_calculated = zero(density)
    for l in 1:num_electron ÷ 2
        density_calculated += 2 .* abs2.(ψ[:,l]) ./ x_step
    end
    if num_electron % 2 == 1
        density_calculated += abs2.(ψ[:,num_electron ÷ 2 + 1]) ./ x_step
    end
    return density_calculated, E, ψ
end

function self_consistent_loop(H::Hamiltonian1D{T},
        num_electron::Int;
        # density_guess = ones(T, H.space.num_samples) .* (num_electron / (H.space.right_boundary - H.space.left_boundary)),
        density_guess = zero(coordinate_sample(H.space)),
        tol = 1e-3, max_iter::Int = 10^4,
        periodic::Bool=true
        ) where T<:Real
    density_old = density_guess
    local density_new
    for _ in 1:max_iter
        density_new,E,_ = kohn_sham_solver(H, num_electron, density_old, periodic=periodic)
        @show E[1]
        (sum(abs2, density_new - density_old) < tol) && break
        density_old = density_new
    end
    return density_new
end

struct NonInteractingHamiltonian1D{RT, FT} <: Hamiltonian1D{RT}
    space::SpaceSample1D{RT}
    external_potential::FT
end

function mat(system::NonInteractingHamiltonian1D{RT, FT},
    density::AbstractVector{RT}; periodic::Bool=true
    ) where {RT<:Real, FT}
return kinetic_energy_1d(system.space,periodic=periodic) +
    external_potential_1d(system.space, system.external_potential)
end

struct InteractingHamiltonian1D{RT, FT} <: Hamiltonian1D{RT}
    space::SpaceSample1D{RT}
    external_potential::FT
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

function one_dim_kronig_penney(a::T, U0::T, space::SpaceSample1D{T}) where T<:Real
    # the Kronig-Penney model is a periodic potential of the form:
    # 1. the period is T = a + b
    # 2. for 0≤x<a, U(x) = 0; for a≤x<a+b, U(x) = U0 
    return NonInteractingHamiltonian1D(space, x -> (x<a ? zero(T) : U0))
end
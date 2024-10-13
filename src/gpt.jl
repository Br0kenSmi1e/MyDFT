using LinearAlgebra

function integral(x, y; axis=1)
    dx = x[2] - x[1]
    return sum(y) * dx
end

function get_nx(num_electron, psi, x)
    I = integral(x, psi.^2)
    normed_psi = psi ./ sqrt(I)
    
    fn = fill(2, num_electron ÷ 2)
    if num_electron % 2 != 0
        push!(fn, 1)
    end
    
    res = zeros(length(x))
    for (ne, psi) in zip(fn, eachcol(normed_psi))
        res .+= ne .* (psi.^2)
    end
    return res
end

function get_exchange(nx, x)
    energy = -3/4 * (3/π)^(1/3) * integral(x, nx.^(4/3))
    potential = -(3/π)^(1/3) * nx.^(1/3)
    return energy, potential
end

function get_hatree(nx, x; eps=1e-1)
    h = x[2] - x[1]
    energy = sum(nx' * nx * h^2 ./ sqrt.((x' .- x).^2 .+ eps) / 2)
    potential = sum(nx' * h ./ sqrt.((x' .- x).^2 .+ eps), dims=2)
    return energy, potential
end

function print_log(i, log)
    println("step: $i energy: $(log["energy"][end]) energy_diff: $(log["energy_diff"][end])")
end

function get_d_d2(h, n_grid)
    # D = -I + Diagonal(ones(n_grid-1), 1)
    D = Tridiagonal(ones(n_grid-1),-2 .* ones(n_grid),ones(n_grid-1))
    D = D / h

    D2 = D * D'
    D2[end, end] = D2[1, 1]
    
    return D, D2
end

# Main code
n_grid = 200
x = range(-5, 5, length=n_grid)
num_electron = 17
max_iter = 1000
energy_tolerance = 1e-5
log = Dict("energy" => [Inf], "energy_diff" => [Inf])

h = x[2] - x[1]
D, D2 = get_d_d2(h, n_grid)

denx = zeros(n_grid)
for i in 0:max_iter
    global denx

    ex_energy, ex_potential = get_exchange(denx, x)
    ha_energy, ha_potential = get_hatree(denx, x)
    # @show size(ex_potential)
    # @show typeof(ex_potential)
    # @show size(ha_potential)
    # @show typeof(ha_potential)
    # @show size(D2)
    # @show typeof(D2)


    H = -D2 ./ 2 + Diagonal(ex_potential .+ ha_potential[1,:] .+ x.^2)
    
    energy, psi = eigen(H)
    
    push!(log["energy"], energy[1])
    energy_diff = energy[1] - log["energy"][end-1]
    push!(log["energy_diff"], energy_diff)
    print_log(i, log)
    
    if abs(energy_diff) < energy_tolerance
        println("converged!")
        break
    end
    
    denx = get_nx(num_electron, psi, x)
end
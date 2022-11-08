include("const.jl")
include("utils.jl")

using LinearAlgebra
# using Plots

function main()
    # if output directory does not exist, create it
    # create_dir("output")
    S = get_s_matrix()
    H = get_kinetic_matrix() + get_parabolic_potential_matrix() + get_gauss_potential_matrix()
    E, c = eigen(H, S)
    println(get_energy_meV(E))
end

function get_kinetic_element(xl, xp, a=A, m=M)
    return -sqrt(a) * exp(-a * (xl - xp)^2 / 2) * (-1 + a * (xl - xp)^2) * sqrt(2) * sqrt(pi) / m / 4
end

function get_parabolic_potential_element(xl, xp, a=A)
    return exp(-a * (xl^2 - 2 * xl * xp + xp^2) / 2) * a^(-1 // 2) * (xl + xp)^2 * sqrt(2) * sqrt(pi) / 8 + exp(-a * (xl^2 - 2 * xl * xp + xp^2) / 2) * sqrt(2) * a^(-3 // 2) * sqrt(pi) / 8
end

function get_gauss_potential_element(xl, xp, a=A, l=L)
    return exp(-a * (l^2 * a * xl^2 - 2 * l^2 * a * xl * xp + l^2 * a * xp^2 + xl^2 + xp^2) / (2 * l^2 * a + 1)) * l * (2 * l^2 * a + 1)^(-1 // 2) * sqrt(pi)
end

function get_s_element(xl, xp, a=A)
    return exp(-a * (xl^2 - 2 * xl * xp + xp^2) / 2) * sqrt(2) * a^(-1 // 2) * sqrt(pi) / 2
end

function get_kinetic_matrix(a=A, m=M)
    kinetic_matrix = zeros(N, N)
    for i in 1:N
        for j in 1:N
            kinetic_matrix[i, j] = get_kinetic_element(i, j, a, m)
        end
    end
    return kinetic_matrix
end

function get_parabolic_potential_matrix(a=A)
    potential_matrix = zeros(N, N)
    for i in 1:N
        for j in 1:N
            potential_matrix[i, j] = get_parabolic_potential_element(i, j, a)
        end
    end
    return potential_matrix
end

function get_gauss_potential_matrix(a=A, l=L)
    potential_matrix = zeros(N, N)
    for i in 1:N
        for j in 1:N
            potential_matrix[i, j] = get_gauss_potential_element(i, j, a, l)
        end
    end
    return potential_matrix
end

function get_s_matrix(a=A)
    s_matrix = zeros(N, N)
    for i in 1:N
        for j in 1:N
            s_matrix[i, j] = get_s_element(i, j, a)
        end
    end
    return s_matrix
end

main()
function get_kinetic_element(xl::Float64, xp::Float64, a::Float64=A, m::Float64=M)::Float64
    return -sqrt(a) * exp(-a * (xl - xp)^2 / 2) * (-1 + a * (xl - xp)^2) * sqrt(2) * sqrt(pi) / m / 4
end

function get_parabolic_potential_element(xl::Float64, xp::Float64, a::Float64=A)::Float64
    return exp(-a * (xl^2 - 2 * xl * xp + xp^2) / 2) * a^(-1 // 2) * (xl + xp)^2 * sqrt(2) * sqrt(pi) / 8 + exp(-a * (xl^2 - 2 * xl * xp + xp^2) / 2) * sqrt(2) * a^(-3 // 2) * sqrt(pi) / 8
end

function get_gauss_potential_element(xl::Float64, xp::Float64, a::Float64=A, l::Float64=L)::Float64
    return exp(-a * (l^2 * a * xl^2 - 2 * l^2 * a * xl * xp + l^2 * a * xp^2 + xl^2 + xp^2) / (2 * l^2 * a + 1)) * l * (2 * l^2 * a + 1)^(-1 // 2) * sqrt(pi)
end

function get_s_element(xl::Float64, xp::Float64, a::Float64=A)::Float64
    return exp(-a * (xl^2 - 2 * xl * xp + xp^2) / 2) * sqrt(2) * a^(-1 // 2) * sqrt(pi) / 2
end

function get_kinetic_matrix(centers::Vector{Float64}, a::Float64=A, m::Float64=M)::Matrix{Float64}
    kinetic_matrix::Matrix{Float64} = zeros(N, N)
    for i::Int64 in 1:N
        for j in 1:N
            kinetic_matrix[i, j] = get_kinetic_element(centers[i], centers[j], a, m)
        end
    end
    return kinetic_matrix
end

function get_parabolic_potential_matrix(centers::Vector{Float64}, a::Float64=A, m::Float64=M, ω::Float64=OMEGA)::Matrix{Float64}
    potential_matrix::Matrix{Float64} = zeros(N, N)
    for i::Int64 in 1:N
        for j::Int64 in 1:N
            potential_matrix[i, j] = get_parabolic_potential_element(centers[i], centers[j], a)
        end
    end
    return m / 2 * ω^2 * potential_matrix
end

function get_gauss_potential_matrix(centers::Vector{Float64}, v0::Float64=V0, a::Float64=A, l::Float64=L)::Matrix{Float64}
    potential_matrix::Matrix{Float64} = zeros(N, N)
    for i::Int64 in 1:N
        for j::Int64 in 1:N
            potential_matrix[i, j] = get_gauss_potential_element(centers[i], centers[j], a, l)
        end
    end
    return v0 * potential_matrix
end

function get_s_matrix(centers::Vector{Float64}, a::Float64=A)::Matrix{Float64}
    s_matrix::Matrix{Float64} = zeros(N, N)
    for i::Int64 in 1:N
        for j::Int64 in 1:N
            s_matrix[i, j] = get_s_element(centers[i], centers[j], a)
        end
    end
    return s_matrix
end
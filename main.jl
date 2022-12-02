include("const.jl")
include("utils.jl")
include("matrices.jl")
include("net.jl")
include("psi.jl")

using LinearAlgebra
using Plots

function get_Energy(dx::Float64, a::Float64, l::Float64, v0::Float64)::Vector{Float64}
    centers::Vector{Float64} = generate_centers(dx)
    S::Matrix{Float64} = get_s_matrix(centers, a)
    K::Matrix{Float64} = get_kinetic_matrix(centers, a)
    V_osc::Matrix{Float64} = get_parabolic_potential_matrix(centers, a)
    V_gauss::Matrix{Float64} = get_gauss_potential_matrix(centers, v0, a, l)
    H::Matrix{Float64} = K + V_osc + V_gauss
    E::Vector{Float64}, c::Matrix{Float64} = eigen(H, S)
    return E
end

function calculate_ψ(dx::Float64, a::Float64, l::Float64, v0::Float64, filename::String)
    centers::Vector{Float64} = generate_centers(dx)
    S::Matrix{Float64} = get_s_matrix(centers, a)
    K::Matrix{Float64} = get_kinetic_matrix(centers, a)
    V_osc::Matrix{Float64} = get_parabolic_potential_matrix(centers, a)
    V_gauss::Matrix{Float64} = get_gauss_potential_matrix(centers, v0, a, l)
    H::Matrix{Float64} = K + V_osc + V_gauss
    E::Vector{Float64}, c::Matrix{Float64} = eigen(H, S)
    println("Energies: ", E * R)
    plot_ψ(centers, c, 1, "$(filename)_1.png")
    plot_ψ(centers, c, 2, "$(filename)_2.png")
    plot_ψ(centers, c, 3, "$(filename)_3.png")
end

function main()
    # if output directory does not exist, create it
    create_dir("output")
    α::Vector{Float64} = (0.1:0.05:2) * A
    dx::Vector{Float64} = (0.1:0.05:2) * Δx
    E_matrix::Matrix{Float64} = zeros(length(α), length(dx))
    for (i, a) in enumerate(α)
        for (j, δx) in enumerate(dx)
            E_matrix[i, j] = get_Energy(δx, a, L, V0)[1]
        end
    end
    heatmap(α / A, dx / Δx, E_matrix * R, xlabel="a", ylabel="dx", title="E_1", color=:thermal)
    savefig("output/E_1.png")


    calculate_ψ(Δx, A, L, V0, "new_psi")

end

main()
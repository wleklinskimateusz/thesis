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

function calculate_ψ(dx::Float64, a::Float64, l::Float64, v0::Float64)
    centers::Vector{Float64} = generate_centers(dx)
    S::Matrix{Float64} = get_s_matrix(centers, a)
    K::Matrix{Float64} = get_kinetic_matrix(centers, a)
    V_osc::Matrix{Float64} = get_parabolic_potential_matrix(centers, a)
    V_gauss::Matrix{Float64} = get_gauss_potential_matrix(centers, v0, a, l)
    H::Matrix{Float64} = K + V_osc + V_gauss
    E, c::Matrix{ComplexF64} = eigen(H, S)
    return c, V_osc + V_gauss
end

function calculate_error(values::Vector{Float64}, teoretical::Vector{Float64})
    error::Vector{Float64} = zeros(length(values))
    for i in 1:length(values)
        error[i] = (values[i] - teoretical[i])^2
    end
    return sum(error)
end

function get_bext_params()
    best_arr = []
    α::Vector{Float64} = (0.1:0.05:3) * A
    dx::Vector{Float64} = (0.1:0.05:3) * Δx
    ε::Float64 = 1e-3
    for a in α
        for δx in dx
            E = get_Energy(δx, a, L, V0)[1:7]
            E_teo = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
            error = calculate_error(E * R, E_teo)
            if error < ε
                push!(best_arr, (δx, a, error))
            end
        end
    end
    return best_arr
end

function plot_probabilities(δx::Float64, α::Float64, states::Int, v0::Float64, title=Nothing, filename=Nothing)
    centers = generate_centers(δx)
    c, V = calculate_ψ(δx, α, L, v0)
    l = maximum(abs.(centers))
    net = generate_net(l)
    V = get_gauss_potential(net, centers, v0, α)
    # V = get_parabolic_potential(net)
    # plot(net, V * R, label="V(x)")
    # savefig("output/V.png")
    p = plot()
    for i in 1:states
        ψ = get_ψ(centers, c, i)
        p = plot!(net, abs2.(ψ), label="ψ$i")
    end
    if (title !== Nothing)
        plot!(title=title)
    end
    if (filename !== Nothing)
        savefig("output/$filename.png")
    end
    return p

end
const eps = 0.0001
function main()
    # if output directory does not exist, create it
    create_dir("output")
    α, δx = 0.55 * A, 0.1 * Δx
    # δx, α, error = get_bext_params()[end]
    α, δx = 1.15 * A, 0.25 * Δx
    anim = Animation()
    for i in 1:100
        v0 = i * 0.01 * V0
        p = plot_probabilities(δx, α, 4, v0, "V0 = $(v0 * R)")
        frame(anim, p)
    end
    gif(anim, "output/PSI.gif", fps=24)
    println("Best parameters: α = $(α / A), δx = $(δx / Δx), error = $error")
    println("Best energy: $(get_Energy(δx, α, L, V0)[1:7] * R)")

end

main()
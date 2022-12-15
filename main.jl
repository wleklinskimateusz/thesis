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
    E::Vector{Float64}, c::Matrix{Float64} = eigen(Hermitian(H), Hermitian(S))
    return E
end

"Returns Matrices: H, V, S in that order"
function generate_matrices(centers::Vector{Float64}, a::Float64, l::Float64, v0::Float64)
    S::Matrix{Float64} = get_s_matrix(centers, a)
    K::Matrix{Float64} = get_kinetic_matrix(centers, a)
    V_osc::Matrix{Float64} = get_parabolic_potential_matrix(centers, a)
    V_gauss::Matrix{Float64} = get_gauss_potential_matrix(centers, v0, a, l)
    V::Matrix{Float64} = V_osc + V_gauss
    return K + V, V, S
end

function calculate_ψ(dx::Float64, a::Float64, l::Float64, v0::Float64)
    centers::Vector{Float64} = generate_centers(dx)
    H, V, S = generate_matrices(centers, a, l, v0)
    E, c::Matrix{Float64} = eigen(Hermitian(H), Hermitian(S))
    return E, c, V
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
    ε::Float64 = 0.1
    for a in α
        for δx in dx
            E = get_Energy(δx, a, L, 0.0)[1:7]
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
    _, c, V = calculate_ψ(δx, α, L, v0)
    l = maximum(abs.(centers))
    net = generate_net(l)
    V = get_gauss_potential(net, centers, v0, α) + get_parabolic_potential(net)
    p = plot()

    for i in 1:states
        ψ = get_ψ(centers, c, i)
        p = plot!(net * L0, 2 * R * abs2.(ψ), label="ψ$i")
    end
    if (title !== Nothing)
        plot!(title=title)
    end
    if (filename !== Nothing)
        savefig("output/$filename.png")
    end
    p = plot!(net * L0, V * R, label="V", xlabel="x [nm]", ylabel="V [meV]")
    return p
end

function animate_v0(δx::Float64, α::Float64, filename::String)
    anim = Animation()
    for i in 1:500
        v0 = i * 0.01 * V0
        p = plot_probabilities(δx, α, 4, v0, "V0 = $(round(v0 * R, digits=2))meV")
        frame(anim, p)
    end
    gif(anim, "output/$(filename).gif", fps=24)
end

function get_next_c_element(c_prev::Vector{Float64}, d::Matrix{Float64}, E0::Vector{Float64}, α::Float64, Hgauss::Matrix{Float64}, l::Int)
    Hz = 0
    n = length(c_prev)
    for k in 1:n
        Hz += c_prev[k] * get_z_element(d, l, k, Hgauss)
    end
    return (1 - E0[l] * α) * c_prev[l] - α * Hz
end

function normalise_eigenvector(c::Vector{Float64}, S::Matrix{Float64})::Vector{Float64}
    norm = c' * S * c
    return c / sqrt(norm)
end

function orthonormalise_eigenvector(c::Vector{Float64}, ci::Vector{Float64})
    r = c' * ci
    return c - r * ci
end

function iterate(δx::Float64, α::Float64, l::Float64, v0::Float64, state::Int, filename::String, prev_c=nothing)
    if (state > 1 && prev_c === nothing)
        error("For excited states you have to pass previous values for c")
    end
    create_dir("output/frames")
    centers = generate_centers(δx)
    E0, d, V = calculate_ψ(δx, α, l, 0.0) # eigenstates for V0 = 0
    H, _, S = generate_matrices(centers, α, l, v0)
    E, c_all = eigen(Hermitian(H), Hermitian(S))
    c::Vector{Float64} = d[:, state]
    println(c' * S * c)
    net = generate_net(l)
    anim = Animation()
    for t in 1:1000000000
        new_c = zeros(length(c))
        for i in 1:N
            new_c[i] = get_next_c_element(c, d, E0, α, H, i)
        end
        if (state > 1)
            for ci in prev_c
                new_c = orthonormalise_eigenvector(new_c, ci)
            end
        end
        c = normalise_eigenvector(new_c, S)
        if (t % 10000000 == 0)
            ψ = zeros(length(net))
            for (i, center) in enumerate(centers)
                ψ += generate_ψ_element(net, c[i], center, α)
            end
            i::Int = trunc(t / 10000000)
            print("$(i)%")
            p = plot(net * L0, abs2.(ψ), xlabel="x [nm]", y="probability [%]", title="t=$t", ylims=(0, 0.002))
            frame(anim, p)
            savefig(p, "output/frames/$i.png")
        end
    end
    gif(anim, "output/$filename.gif")
    return c
end

function save_c_to_file(c::Vector{Float64}, filename::String)
    create_dir("output/vectors")
    f = open("output/vectors/$filename.txt", "w")
    for (i, element) in enumerate(c)
        if (i == length(c))
            write(f, "$element")
        else
            write(f, "$element ")
        end
    end
    close(f)
end

function read_from_file(filename::String)::Vector{Float64}
    f = open("output/vectors/$filename.txt")
    file = read(f, String)
    elements::Vector{Float64} = [parse(Float64, element) for element in split(file, " ")]
    close(f)
    return elements
end

const eps = 0.0001
function main()
    # if output directory does not exist, create it
    create_dir("output")
    # α, δx = 0.55 * A, 0.1 * Δx
    # δx, α, error = @time get_bext_params()[end]
    α, δx = 2.55 * A, 0.25 * Δx
    # @time animate_v0(δx, α, "Psi2")
    c1 = @time iterate(δx, α, L, V0, 1, "iteration")
    save_c_to_file(c1, "c1")
    c1 = read_from_file("c1")
    # c2 = @time iterate(δx, α, L, V0, 2, "Psi2", [c1])
    # save_c_to_file(c2, "c2")


end

main()
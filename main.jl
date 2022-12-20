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
function generate_matrices(centers::Vector{Float64}, a::Float64, l::Float64, v0::Float64, include_gauss::Bool=true)
    S::Matrix{Float64} = get_s_matrix(centers, a)
    K::Matrix{Float64} = get_kinetic_matrix(centers, a)
    V_osc::Matrix{Float64} = get_parabolic_potential_matrix(centers, a)
    V_gauss::Matrix{Float64} = get_gauss_potential_matrix(centers, v0, a, l)
    V::Matrix{Float64} = V_osc
    if include_gauss
        V += V_gauss
    end
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
    l = 150 / L0
    net = generate_net(l)
    V = get_gauss_potential(net, centers, v0, α) + get_parabolic_potential(net)
    p = plot()

    for i in 1:states
        ψ = get_ψ(net, centers, c, i)
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

function get_states(δx::Float64, α::Float64, states::Int, v0::Float64)
    centers = generate_centers(δx)
    _, c, V = calculate_ψ(δx, α, L, v0)
    l = 150 / L0
    net = generate_net(l)
    V = get_parabolic_potential(net)
    states::Matrix{Float64} = zeros(length(net), states)
    for i in 1:states
        ψ = get_ψ(centers, c, i)
        states[:, i] = ψ
    end
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

function get_next_c_element(centers::Vector{Float64}, c_prev::Vector{Float64}, d::Matrix{Float64}, E0::Vector{Float64}, α::Float64, l::Int, v0::Float64)
    Hz = 0
    iter_param = 10
    for k in 1:N
        Hz += c_prev[k] * get_z_element(centers, d, l, k, α, v0)
    end
    return (1 - E0[l] * iter_param) * c_prev[l] - iter_param * Hz
end

function normalise_eigenvector(c::Vector{Float64})::Vector{Float64}
    norm = c' * c
    return c / sqrt(norm)
end

function orthonormalise_eigenvector(c::Vector{Float64}, ci::Vector{Float64})
    r = c' * ci
    return c - r * ci
end

function getEnergy(centers::Vector{Float64}, d::Matrix{Float64}, c::Vector{Float64}, E0::Vector{Float64}, α::Float64, v0::Float64)
    output::Float64 = 0
    for k in 1:N
        output += c[k]^2 * E0[k]
        for l in 1:N
            Z = get_z_element(centers, d, l, k, α, v0)
            output += c[k] * c[l] * Z
        end
    end
    return output
end

function iterate(δx::Float64, α::Float64, l::Float64, v0::Float64, state::Int, filename::String, prev_c=nothing, init_c=nothing, t0=1, t_num=100000000)
    if (state > 1 && prev_c === nothing)
        error("For excited states you have to pass previous values for c")
    end
    create_dir("output/frames")
    centers = generate_centers(δx)
    E0, d, V = calculate_ψ(δx, α, l, 0.0) # eigenstates for V0 = 0
    _, ddiag, _ = calculate_ψ(δx, α, l, v0) # eigenstates for V0 = v0
    H, _, S = generate_matrices(centers, α, l, v0, false)
    # E, c_all = eigen(Hermitian(H), Hermitian(S))
    c::Vector{Float64} = zeros(N)
    if (init_c === nothing)
        c[1] = 1
    else
        c = init_c
    end
    println(c' * c)
    x::Vector{Float64} = (-200:200) / L0
    base = zeros((length(x), N))
    for (i, x_i) in enumerate(x)
        for j in 1:N
            output = 0
            for k in 1:N
                xk = (k - N / 2) * δx
                output += H[k, j] * exp(-α * (x_i - xk)^2)
            end
            base[i, j] = output
        end
    end
    plot(x * L0, base)
    savefig("output/state_1.png")
    anim = Animation()
    E::Vector{Float64} = []
    for t in t0:t0+t_num
        new_c = zeros(N)
        for i in 1:N
            new_c[i] = get_next_c_element(centers, c, d, E0, α, i, v0)
        end
        if (state > 1)
            for ci in prev_c
                new_c = orthonormalise_eigenvector(new_c, ci)
            end
        end
        c = normalise_eigenvector(new_c)
        freq::Int = trunc(t_num / 100)
        if (t % freq == 0)
            ψ = zeros(length(x))

            for i in 1:N
                ψ += c[i] * get_ψ(x, centers, d, i)
            end

            i::Int = trunc(t / freq)
            print("$i%", "\r\r\r")
            p = plot(x * L0, abs2.(ψ), xlabel="x [nm]", y="probability [%]", title="t=$t")
            push!(E, getEnergy(centers, d, c, E0, α, v0))
            frame(anim, p)
            savefig(p, "output/frames/$i.png")
            save_c_to_file(c, "temp/c")
        end
    end
    gif(anim, "output/$filename.gif")
    ψ = zeros(length(x))
    diagState = state == 2 ? 3 : 1
    ψdiag = generate_ψ(x, centers, ddiag[:, diagState])
    for i in 1:N
        ψ += c[i] * generate_ψ(x, centers, d[:, i])
    end

    plot(x * L0, abs2.(ψ), xlabel="x [nm]", y="probability [%]", label="Imaginary time evolution")
    plot!(x * L0, abs2.(ψdiag), xlabel="x [nm]", y="probability [%]", label="Diagonalisation")
    savefig("output/comparing$(state).png")

    save_c_to_file(E, "E")
    p = plot(E * R)
    savefig(p, "output/Energy.png")
    return c
end

function save_c_to_file(c::Vector{Float64}, filename::String)
    create_dir("output/vectors")
    create_dir("output/vectors/temp")
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
    # print(get_Energy(δx, α, L, 0.0) * R)
    # α, δx = 2.55 * A, 0.25 * Δx
    α = 2.55 * A
    δx = 0.25 * Δx
    println("δx = $(δx * L0), α = $α")
    # @time animate_v0(δx, α, "Psi2")
    # c1 = read_from_file("c1")
    c1 = @time iterate(δx, α, L, V0, 1, "iteration", nothing, nothing, 1, 10000)
    save_c_to_file(c1, "c1")
    c2 = @time iterate(δx, α, L, V0, 2, "iteration2", [c1], nothing, 1, 10000)
    save_c_to_file(c2, "c2")
    # c3 = @time iterate(δx, α, L, V0, 3, "iteration3", [c1, c2], nothing, 1, 10000)
    # save_c_to_file(c3, "c3")
    # c2 = @time iterate(δx, α, L, V0, 2, "Psi2", [c1])
    # save_c_to_file(c2, "c2")
    # c3 = @time iterate(δx, α, L, V0, 3, "Psi3", [c1, c2])
    # save_c_to_file(c3, "c3")
    # c4 = @time iterate(δx, α, L, V0, 4, "Psi4", [c1, c2, c3])
    # save_c_to_file(c4, "c4")


end

main()
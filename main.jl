include("const.jl")
include("utils.jl")
include("matrices.jl")
include("net.jl")
include("psi.jl")



using LinearAlgebra
using Plots;
font = Plots.font("Helvetica", 14)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

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
    E::Vector{Float64}, c::Matrix{Float64} = eigen(Hermitian(H), Hermitian(S))
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
    α::Vector{Float64} = (0.5:0.005:3) * A
    dx::Vector{Float64} = (0.05:0.001:0.4) * Δx
    ε::Matrix{Float64} = zeros(length(α), length(dx))
    for (i, a) in enumerate(α)
        for (j, δx) in enumerate(dx)
            E = get_Energy(δx, a, L, 0.0)
            E_teo = [get_E_teo(i) for i in 1:9]
            error = calculate_error(E * R, E_teo)
            ε[i, j] = error
        end
    end
    heatmap(dx * L0, α / A, log.(10, ε), ylabel="a / mω", xlabel="Δx [nm]", legend=true)
    savefig("output/R2.png")
    # save 5 best parameters
    for i in 1:10
        min = findmin(ε)
        ε[min[2][1], min[2][2]] = Inf
        push!(best_arr, (α[min[2][1]] / A, dx[min[2][2]] * L0))
    end

    return best_arr
end

function plot_gauss_base(α::Float64, dx::Float64, l::Float64, filename::String)
    centers::Vector{Float64} = generate_centers(dx)
    net = generate_net(150 / L0)
    output::Vector{Float64} = zeros(length(net))
    plot(xlabel="x [nm]", ylabel="A.U.", title="Gauss base")
    for center in centers
        plot!(net * L0, exp.(-α * (net .- center) .^ 2), label="")
    end
    savefig("output/$filename.png")
end

function plot_probabilities(δx::Float64, α::Float64, states::Int, v0::Float64, title=Nothing, filename=Nothing, plot_energy=false)
    centers = generate_centers(δx)
    E, c, V = calculate_ψ(δx, α, L, v0)
    if (plot_energy)
        scatter(1:states, E[1:states] * R, xlabel="N", ylabel="E [meV]", label=false)
        savefig("output/E_$filename.png")
    end

    l = 200 / L0
    net = generate_net(l)
    V = get_gauss_potential(net, v0, α) + get_parabolic_potential(net)
    p = plot()

    for i in 1:states
        ψ = get_ψ(net, centers, c, i)
        p = plot!(net * L0, abs2.(ψ), label="N=$(i-1)", xlabel="x [nm]", ylim=(0, 0.001))
    end
    if (title !== Nothing)
        plot!(title=title)
    end
    if (filename !== Nothing)
        savefig("output/$filename.png")
    end
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
        ψ = get_ψ(net, centers, c, i)
        states[:, i] = ψ
    end
end

function animate_v0(δx::Float64, α::Float64, filename::String)
    to_save = [0.03, 0.51, 0.99, 2.01, 3.0, 3.99, 5.01, 9.00, 15.0]
    anim = Animation()
    for i in 1:500
        v0 = i * 0.01 * V0
        p = plot_probabilities(δx, α, 3, v0, "V0 = $(round(v0 * R, digits=2))meV")
        if (round(v0 * R, digits=2) in to_save)
            savefig("output/anim/frames/$(round(v0 * R, digits=2)).png")
        end
        frame(anim, p)
    end
    gif(anim, "output/$(filename).gif", fps=24)
end

function get_next_c_element(centers::Vector{Float64}, c_prev::Vector{Float64}, d::Matrix{Float64}, E0::Vector{Float64}, α::Float64, l::Int, v0::Float64)
    Hz = 0
    iter_param = 15
    for k in 1:N
        Hz += c_prev[k] * get_z_element(centers, d, l, k, α, v0)
    end
    return (1 - E0[l] * iter_param) * c_prev[l] - iter_param * Hz
end

function normalise_eigenvector(c::Vector{Float64})::Vector{Float64}
    norm = c' * c
    return c / sqrt(norm)
end

function orthogonalise(c::Vector{Float64}, ci::Vector{Float64})
    r = c' * ci
    return c - r * ci
end

function calculate_energy(centers::Vector{Float64}, d::Matrix{Float64}, c::Vector{Float64}, E0::Vector{Float64}, α::Float64, v0::Float64)
    E::Float64 = 0
    for k in 1:N
        E += c[k]^2 * E0[k]
        for l in 1:N
            Z = get_z_element(centers, d, l, k, α, v0)
            E += c[k] * c[l] * Z
        end
    end
    return E
end


function iterate(δx::Float64, α::Float64, l::Float64, v0::Float64, state::Int, filename::String, prev_c=nothing, init_c=nothing, t0=1, t_num=100000000)
    if (state > 1 && prev_c === nothing)
        error("For excited states you have to pass previous values for c")
        if (length(prev_c) != N - 1)
            error("Length of previous c is not enough to calculate this state")
        end
    end
    create_dir("output/frames")
    centers = generate_centers(δx)
    E0, d, V = calculate_ψ(δx, α, l, 0.0) # eigenstates for V0 = 0
    Ediag, ddiag, _ = calculate_ψ(δx, α, l, v0) # eigenstates for V0 = v0
    H, _, S = generate_matrices(centers, α, l, v0, false)
    # E, c_all = eigen(Hermitian(H), Hermitian(S))
    c::Vector{Float64} = zeros(N)
    c[state] = 1
    # if (state % 2 == 1)
    #     c[1] = 1
    # else
    #     c[1] = 0.5
    #     c[2] = 0.5
    # end
    c = normalise_eigenvector(c)
    x::Vector{Float64} = (-200:200) / L0
    anim = Animation()
    E::Vector{Float64} = []
    freq::Int = trunc(t_num / 100)
    for t in t0:t0+t_num
        if (t % freq == 0 || t == t0)
            ψ = zeros(length(x))
            for i in 1:N
                ψ += c[i] * generate_ψ(x, centers, d[:, i])
            end

            i::Int = trunc(t / freq)
            print("$i%", "\r\r\r")
            p = plot(x * L0, abs2.(ψ), xlabel="x [nm]", ylabel="probability density", label="ψ$(state)", title="t=$t", ylims=(0, 0.0012))
            push!(E, calculate_energy(centers, d, c, E0, α, v0))
            frame(anim, p)
            savefig(p, "output/frames/$i.png")
            save_c_to_file(c, "temp/c")
        end
        new_c = zeros(N)
        for i in 1:N
            new_c[i] = get_next_c_element(centers, c, d, E0, α, i, v0)
        end
        if (state > 1)
            for ci in prev_c
                new_c = orthogonalise(new_c, ci)
            end
        end
        c = normalise_eigenvector(new_c)
    end
    gif(anim, "output/$filename.gif")
    ψ = zeros(length(x))
    ψdiag = get_ψ(x, centers, ddiag, state)
    for i in 1:N
        ψ += c[i] * generate_ψ(x, centers, d[:, i])
    end
    ψ = normalise(x, ψ)

    plot(x * L0, abs2.(ψ), xlabel="x [nm]", y="probability [%]", label="Imaginary time evolution", ylims=(0, 0.0015))
    plot!(x * L0, abs2.(ψdiag), xlabel="x [nm]", y="probability [%]", label="Diagonalisation", ylims=(0, 0.0015))
    savefig("output/comparing$(state).png")

    save_c_to_file(E, "E")
    t::Vector{Int} = round.(Int, LinRange(t0, t0 + t_num, length(E)))

    plot(t, Ediag[state] * R * ones(length(t)), xlabel="iteration step", ylabel="Energy [meV]", label="Diagonalisation", formatter=:plain, linestyle=:dash)
    plot!(t, E * R, xlabel="iteration step", ylabel="Energy [meV]", label="Imaginary time evolution", formatter=:plain)
    savefig("output/Energy_state$(state).png")
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

function get_E_teo(state::Int)::Float64
    return (state - 1 / 2) * OMEGA * R
end

function compare_energies(E::Vector{Float64})
    E_teo = [get_E_teo(state) for state in 1:length(E)]
    states = 1:length(E)
    scatter(states, E, xlabel="stany", ylabel="E [meV]", label="wartości numeryczne")
    scatter!(states, E_teo, label="teoretyczne wartości")

    savefig("output/E_teo_vs_E.png")
end

function plot_potential()
    x::Vector{Float64} = (-200:200) / L0
    Vh = get_parabolic_potential(x)
    Vz = get_gauss_potential(x, V0, L)
    V = Vh + Vz
    plot(x * L0, V * R, xlabel="x [nm]", ylabel="V [meV]", label="V(x)", ylim=(0, 15))
    savefig("output/V.png")
    plot(x * L0, Vh * R, xlabel="x [nm]", ylabel="V [meV]", label="Vh(x)", ylim=(0, 15))
    savefig("output/Vh.png")
    plot(x * L0, Vz * R, xlabel="x [nm]", ylabel="V [meV]", label="Vz(x)", ylim=(0, 15))
    savefig("output/Vz.png")
end

function plot_energy(δx::Float64, a::Float64, L::Float64, v0::Float64)
    E, _, _ = calculate_ψ(δx, a, L, v0)
    E0, _, _ = calculate_ψ(δx, a, L, 0.0)
    states = 1:length(E)
    scatter(states, E * R, xlabel="stany", ylabel="E [meV]", label="Z zaburzeniem")
    scatter!(states, E0 * R, label="Bez zaburzenia")
    savefig("output/Energy_compare.png")
end

const eps = 0.0001
function main()
    # if output directory does not exist, create it
    create_dir("output")
    create_dir("output/anim")
    create_dir("output/anim/frames")
    plot_potential()
    # α, δx = 0.55 * A, 0.1 * Δx
    # params = @time get_bext_params()
    # println(params)
    # α, δx = 2.55 * A, 0.25 * Δx
    α = 2.55 * A
    δx = 0.25 * Δx
    # α = A / 4
    # δx = 20 / L0
    α = A * 1.5
    δx = 25.5966 / L0
    # α = A * 0.52
    # δx = 7.89 / L0
    println("δx = $(δx * L0), α = $α")
    # print(get_Energy(δx, α, L, 0.0) * R)
    # compare_energies(get_Energy(δx, α, L, 0.0) * R)

    # plot_gauss_base(α, δx, L, "gauss_base")
    # println("δx = $(δx * L0), α = $α")
    # @time plot_energy(δx, α, L, V0)
    # # @time plot_probabilities(δx, α, 3, 0.0, "Oscylator Harmoniczny", "harmonic_oscilator")
    # # @time plot_probabilities(δx, α, 3, V0, "V0=3meV", "with_potential")
    # @time animate_v0(δx, α, "Psi2")
    c1 = read_from_file("c1")
    # c1 = @time iterate(δx, α, L, V0, 1, "iteration", nothing, nothing, 1, 13000)
    # save_c_to_file(c1, "c1")
    # c1 = read_from_file("c1")
    # c2 = @time iterate(δx, α, L, V0, 2, "iteration2", [c1], nothing, 1, 13000)
    # save_c_to_file(c2, "c2")
    c2 = read_from_file("c2")
    # c3 = @time iterate(δx, α, L, V0, 3, "iteration3", [c1, c2], nothing, 1, 13000)
    # save_c_to_file(c3, "c3")
    c3 = read_from_file("c3")
    c4 = @time iterate(δx, α, L, V0, 4, "iteration4", [c1, c2, c3], nothing, 1, 13000)
    save_c_to_file(c4, "c4")
    # c4 = read_from_file("c4")
    c5 = @time iterate(δx, α, L, V0, 5, "iteration5", [c1, c2, c3, c4], nothing, 1, 13000)
    save_c_to_file(c5, "c5")
    # c5 = read_from_file("c5")
    c6 = @time iterate(δx, α, L, V0, 6, "iteration6", [c1, c2, c3, c4, c5], nothing, 1, 13000)
    save_c_to_file(c6, "c6")
    # c6 = read_from_file("c6")
    c7 = @time iterate(δx, α, L, V0, 7, "iteration7", [c1, c2, c3, c4, c5, c6], nothing, 1, 13000)
    save_c_to_file(c7, "c7")
    # c7 = read_from_file("c7")
    c8 = @time iterate(δx, α, L, V0, 8, "iteration8", [c1, c2, c3, c4, c5, c6, c7], nothing, 1, 13000)
    save_c_to_file(c8, "c8")
    # c8 = read_from_file("c8")
    c9 = @time iterate(δx, α, L, V0, 9, "iteration9", [c1, c2, c3, c4, c5, c6, c7, c8], nothing, 1, 13000)
    save_c_to_file(c9, "c9")
    # # c2 = @time iterate(δx, α, L, V0, 2, "Psi2", [c1])
    # # save_c_to_file(c2, "c2")
    # # c3 = @time iterate(δx, α, L, V0, 3, "Psi3", [c1, c2])
    # # save_c_to_file(c3, "c3")
    # # c4 = @time iterate(δx, α, L, V0, 4, "Psi4", [c1, c2, c3])
    # # save_c_to_file(c4, "c4")
    # c1 = read_from_file("c1")
    # c2 = read_from_file("c2")
    # c3 = read_from_file("c3")
    # centers = generate_centers(δx)
    # E0, d, V = calculate_ψ(δx, α, L, 0.0)
    # x::Vector{Float64} = (-200:200) / L0
    # plot()
    # for c in [c1, c2, c3]
    #     ψ = zeros(Float64, length(x))
    #     for i in 1:N
    #         ψ += c[i] * get_ψ(x, centers, d, i)
    #     end
    #     plot!(x, abs2.(ψ))
    # end
    # savefig("output/abs2psi.png")


end

main()
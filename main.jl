include("const.jl")
include("utils.jl")
include("matrices.jl")
include("net.jl")
include("psi.jl")

using LinearAlgebra
using Plots

function main()
    # if output directory does not exist, create it
    create_dir("output")
    centers::Vector{Float64} = generate_centers()
    S::Matrix{Float64} = get_s_matrix(centers)
    K::Matrix{Float64} = get_kinetic_matrix(centers)
    V_osc::Matrix{Float64} = get_parabolic_potential_matrix(centers)
    V_gauss::Matrix{Float64} = get_gauss_potential_matrix(centers)

    H = K + V_osc + V_gauss
    E, c = eigen(H, S)


    println(get_energy_meV(E))
    save_energy_to_file(E, "energy")
    plot_ψ(centers, c, 1, "psi1")
    plot_ψ(centers, c, 2, "psi2")
    plot_ψ(centers, c, 3, "psi3")
    plot_ψ(centers, c, 4, "psi4")

end

main()
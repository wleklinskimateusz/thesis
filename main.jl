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
    centers = generate_centers()
    S = get_s_matrix(centers)
    K = get_kinetic_matrix(centers)
    V_osc = get_parabolic_potential_matrix(centers)
    V_gauss = get_gauss_potential_matrix(centers)

    H = K + V_osc + V_gauss
    E, c = eigen(H, S)


    println(get_energy_meV(E))
    save_energy_to_file(E, "energy")
    plot_ψ(centers, c, 1, "psi1")

end

main()
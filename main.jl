include("const.jl")
include("utils.jl")
include("matrices.jl")
include("net.jl")

using LinearAlgebra
# using Plots

function main()
    # if output directory does not exist, create it
    create_dir("output")
    net = generate_net()
    S = get_s_matrix(net)
    K = get_kinetic_matrix(net)
    V_osc = get_parabolic_potential_matrix(net)
    V_gauss = get_gauss_potential_matrix(net)

    H = K + V_osc + V_gauss
    E, c = eigen(H, S)

    println(get_energy_meV(E))
    save_energy_to_file(E, "energy")

end

main()
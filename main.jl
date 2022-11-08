include("const.jl")
include("utils.jl")
include("matrices.jl")
include("net.jl")

using LinearAlgebra
# using Plots

function main()
    # if output directory does not exist, create it
    # create_dir("output")
    net = generate_net()
    S = get_s_matrix(net)
    H = get_kinetic_matrix(net) + get_parabolic_potential_matrix(net) + get_gauss_potential_matrix(net)
    E, c = eigen(H, S)
    println(get_energy_meV(E))
end

main()
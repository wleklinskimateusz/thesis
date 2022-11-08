include("const.jl")
include("utils.jl")

using LinearAlgebra
# using Plots

function main()
    # if output directory does not exist, create it
    # create_dir("output")

    S = get_s_matrix()
    H = get_kinetic_matrix() + get_parabolic_potential_matrix() + get_gauss_potential_matrix()
    E, c = eigen(H, S)
    println(E)
end

main()
function create_dir(name::String)::Nothing
    if !isdir(name)
        mkdir(name)
    end
    return nothing
end

function get_energy_meV(energy::Float64)::Float64
    return energy * R
end

function get_energy_atomic(energy::Float64)::Float64
    return energy / (R)
end

function save_energy_to_file(energy::Vector{Float64}, filename::String)::Nothing

    open("output/$filename.txt", "w") do f
        for e::Float64 in energy
            write(f, "$(get_energy_meV(e)) ")
        end
    end
end

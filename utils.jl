function create_dir(name)
    if !isdir(name)
        mkdir(name)
    end
end

function get_energy_meV(energy)
    return energy * 2 / R
end

function get_energy_atomic(energy)
    return energy / (2 * R)
end

function save_energy_to_file(energy, filename)

    open("output/$filename.txt", "w") do f
        for e in energy
            write(f, "$(get_energy_meV(e)) ")
        end
    end
end

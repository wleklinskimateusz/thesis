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

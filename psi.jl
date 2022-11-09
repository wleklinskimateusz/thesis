function generate_ψ_element(net, c, xn, a=A)
    return c .* exp.(-a * (net .- xn) .^ 2)
end

function generate_ψ(net, centers, c, a=A)
    output = zeros(length(net))
    for (center, ci) in zip(centers, c)
        output += generate_ψ_element(net, ci, center, a)
    end
    return output
end

function normalise(ψ)
    return ψ ./ sqrt(sum(ψ .^ 2))
end

function plot_ψ(centers, c, i, filename=Nothing)
    net = generate_net()
    ψ = generate_ψ(net, centers, c, i)
    ψ = normalise(ψ)
    plot(net, ψ, label="ψ$i")
    if filename != Nothing
        savefig("output/$filename.png")
    else
        savefig("output/psi$i.png")
    end
end
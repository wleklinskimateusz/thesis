function generate_ψ_element(net::Vector{Float64}, c::Float64, xn::Float64, a::Float64=A)::Vector{Float64}
    return c .* exp.(-a * (net .- xn) .^ 2)
end

function generate_ψ(net::Vector{Float64}, centers::Vector{Float64}, c::Vector{Float64}, a::Float64=A)::Vector{Float64}
    output::Vector{Float64} = zeros(length(net))
    for (center::Float64, ci::Float64) in zip(centers, c)
        output += generate_ψ_element(net, ci, center, a)
    end
    return output
end

function normalise(net::Vector{Float64}, ψ::Vector{Float64})::Vector{Float64}
    dx = net[2] - net[1]
    return ψ ./ sqrt(sum(ψ .^ 2) * dx)
end

function get_ψ(x::Vector{Float64}, centers::Vector{Float64}, c::Matrix{Float64}, i::Int64)::Vector{Float64}
    return normalise(x, generate_ψ(x, centers, c[:, i]))
end

function plot_ψ(centers::Vector{Float64}, c::Matrix{Float64}, i::Int64, filename=Nothing)
    # find max value of centers
    l = maximum(abs.(centers))
    net = generate_net(l)

    ψ = get_ψ(net, centers, c, i)
    plot(net, abs2.(ψ), label="ψ$i")
    if filename != Nothing
        savefig("output/$filename.png")
    else
        savefig("output/psi$i.png")
    end
    return nothing
end
function generate_centers(dx::Float64, n::Int64=N ÷ 2)::Vector{Float64}
    # generate centers with step dx and n centers on each side with one at 0
    start::Float64 = -n * dx
    stop::Float64 = n * dx
    return start:dx:stop
end

function generate_net(l::Float64)::Vector{Float64}
    start = -l / 2
    stop = l / 2
    step = l / 500
    return start:step:stop
end
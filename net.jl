"generate centers with step dx and n centers on each side with one at 0"
function generate_centers(dx::Float64, n::Int64=N ÷ 2)::Vector{Float64}
    start::Float64 = -n * dx
    stop::Float64 = n * dx
    output::Vector{Float64} = start:dx:stop
    return output
end

function generate_net(l::Float64)::Vector{Float64}
    start = -l
    stop = l
    step = l / 100
    return start:step:stop
end
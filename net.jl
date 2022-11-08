function generate_net(l=L, n=N)
    start = -l / 2
    stop = l / 2
    step = l / (n - 1)
    return start:step:stop
end
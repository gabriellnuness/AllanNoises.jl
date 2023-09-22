using AllanNoises
using AllanDeviations
using Test
using PyPlot
using Random
using Printf



@testset "Random Walk" begin
    fs = 100.0
    time = 5*60*60
    t = range(start=0, stop=time, step=1/fs)
    n = length(t)

    Random.seed!(3)
    y = @. 2*randn(n) + cumsum(2e-4*randn(n))
    adev = allandev(y, fs, frequency=true)

    noises = allan_fit(adev.tau, adev.deviation)
    close("all");figure(); 
    yscale("log"); xscale("log")
    plot(adev.tau, adev.deviation);
    plot(adev.tau, noises[1]);
    noises_string = @sprintf "Q = %.2e\narw = %.2e\ndrift = %.2e\nrrw = %.2e\nrr = %.2e" noises[2] noises[3] noises[4] noises[5] noises[6] 
    legend(["data",noises_string])
end
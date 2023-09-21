using AllanFit
using AllanDeviations
using Test
using PyPlot


y = randn(1000)
adev = allandev(y,1.0)

loglog(adev.tau, adev.deviation)


@testset "Random Walk" begin


    @test 1 == 1 # dumb test
    noises = allan_fit(adev.tau, adev.deviation)


end
using AllanNoises
using AllanDeviations
using Test
using PyPlot
using Random
using Printf



@testset "Default signal" begin
    
    Random.seed!(3)
    # creating signal
    fs = 100.0
    time = 5*60*60
    t = range(start=0, stop=time, step=1/fs)
    n = length(t)
    y = @. 2*randn(n) + cumsum(2e-4*randn(n))   # creating a simulated noisy signal
    adev = allandev(y, fs, frequency=true)

    noises = allan_fit(adev.tau, adev.deviation)
    noises_simp = allan_fit_simple(adev.tau, adev.deviation)
    

    
    @test noises.arw    ≈ noises_simp.arw       atol=1e-3
    @test noises.bias   ≈ noises_simp.bias      atol=3e-2   # bias drift is always distant with this method
    @test noises.rrw    ≈ noises_simp.rrw       atol=1e-3



    # Plot fitting
    close("all");figure(); 
    yscale("log"); xscale("log")
    plot(adev.tau, adev.deviation);
    plot(adev.tau, noises[1]);
    noises_string = @sprintf "Q = %.2e\narw = %.2e\ndrift = %.2e\nrrw = %.2e\nrr = %.2e" noises[2] noises[3] noises[4] noises[5] noises[6] 
    legend(["data",noises_string])

    # Plot simple fitting
    lineN = noises_simp.arw ./ sqrt.(adev.tau)
    lineB = noises_simp.bias * 0.664 * ones(size(adev.tau));
    lineK = noises_simp.rrw .* sqrt.(adev.tau/3);
    figure()
    loglog(adev.tau, adev.deviation, color="tab:gray")
    loglog(adev.tau, lineN,"--",color="tab:blue",linewidth=.8)
    loglog(1,noises_simp.arw,"o",color="tab:blue")
    loglog(adev.tau, lineB,"--",color="tab:purple",linewidth=.8)
    loglog(adev.tau[noises_simp.ind_bias], noises_simp.bias*0.664,"o",color="tab:purple")
    loglog(adev.tau, lineK,"--",color="tab:green",linewidth=.8)
    loglog(3,noises_simp.rrw,"o",color="tab:green")
    noises_string = @sprintf "arw = %.2e\ndrift = %.2e\nrrw = %.2e" noises_simp.arw noises_simp.bias noises_simp.rrw
    legend(["data",noises_string])

end
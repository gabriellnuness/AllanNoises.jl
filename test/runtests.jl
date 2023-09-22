using AllanNoises
using AllanDeviations
using Test
using PyPlot
using Random
using Printf

custom_plot_colors =   ["#4063d8"
                        "#9558b2"
                        "#389826"
                        "#cb3c33"
                        "#17becf"
                        "#e377c2"
                        "#8c564b"
                        "#7f7f7f"
                        "#bcbd22"]

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
    close("all");
    figure(); 
    yscale("log"); xscale("log")
    plot(adev.tau, adev.deviation);
    plot(adev.tau, noises[1], color=custom_plot_colors[4]);
        noises_string = @sprintf "Q = %.2e\narw = %.2e\ndrift = %.2e\nrrw = %.2e\nrr = %.2e" noises[2] noises[3] noises[4] noises[5] noises[6] 
        legend(["data",noises_string])
        xlabel("Correlation time [s]")
        ylabel("Signal amplitude [a.u.]")

    # Plot simple fitting
    line_arw = noises_simp.arw ./ sqrt.(adev.tau)
    line_bias = noises_simp.bias * 0.664 * ones(size(adev.tau));
    line_rrw = noises_simp.rrw .* sqrt.(adev.tau/3);
    figure()
    loglog(adev.tau, adev.deviation)
    loglog(adev.tau, line_arw,"--", color=custom_plot_colors[2],linewidth=.8, label="_nolegend_")
    loglog(1, noises_simp.arw,"o",markersize=3,color=custom_plot_colors[2])
    loglog(adev.tau, line_bias,"--",color=custom_plot_colors[3],linewidth=.8, label="_nolegend_")
    loglog(adev.tau[noises_simp.ind_bias], noises_simp.bias*0.664,"o", markersize=3, color=custom_plot_colors[3])
    loglog(adev.tau, line_rrw,"--",color=custom_plot_colors[4],linewidth=.8, label="_nolegend_")
    loglog(3, noises_simp.rrw,"o",markersize=3,color=custom_plot_colors[4])
        str_arw = @sprintf "arw = %.2e" noises_simp.arw
        str_bias = @sprintf "drift = %.2e" noises_simp.bias
        str_rrw = @sprintf "rrw = %.2e" noises_simp.rrw
        legend(["data",str_arw,str_bias,str_rrw])
        xlabel("Correlation time [s]")
        ylabel("Signal amplitude [a.u.]")


end
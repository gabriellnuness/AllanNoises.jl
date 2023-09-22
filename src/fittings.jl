
"""
    Allan noises fitting

 Calculate the sensor noise parameters/coefficients from the 
 Allan Variance data using Weighted Linear Regression

Fitted Allan deviation noise values:

 C = avar_fit
 unit = measurement unit in SI

 Q     --> unit * s
 arw   --> unit * s^(1/2)
 bias  --> unit * s^(0)
 rrw   --> unit * s^(-1/2) 
 rr    --> unit * s^(-1)

 ref: 2017_Jerath-Bridging the gap between sensor noise modeling and sensor characterization (https://doi.org/10.1016/j.measurement.2017.09.012 )
"""
function allan_fit(taus::Array{Float64}, adev::Array{Float64})
    
    ## Weighted Least-Squares Fitting
    avar = adev.^2

    # The weight is approximately to the order of magnitude
    weight = 1 ./avar
    W = diagm(weight)
    T = [taus.^(-2) taus.^(-1) taus.^(0) taus.^(1) taus.^(2)]

    A = T'*W*T
    b = T'*W*avar
    prob = LinearProblem(A, b)
    avar_fit = solve(prob)

    # Check non trustable values from fit
    for i in eachindex(avar_fit) 
        if avar_fit[i] < 0
            avar_fit[i] = 0
            @warn("Noise $i value did not fit correctly")
        end
    end

    Q    = (avar_fit[1]/3)^(0.5)    # C = 3*Q^2 Q^2 = (q^2)/12, where q = quantum step size
    arw  = avar_fit[2]^(0.5)                        # C = N^2
    bias = (avar_fit[3]^(0.5)) / sqrt(2*log(2)/pi)  # C = (0.6643*B)^2
    rrw  = (3*avar_fit[4])^(0.5)                    # C = K^2 /3
    rr   = (2*avar_fit[5])^(0.5)                    # C = R^2 / 2

    avar_fit = sqrt.(avar_fit'*T')
    avar_fit = avar_fit[1,:]

    return (fit=avar_fit, q=Q, arw=arw, bias=bias, rrw=rrw, rr=rr)
end






"""
    Fittins one slope at a time

Fits ARW, RRW, and Bias drift from Allan deviation plot
arw,  bias,   rrw       = Fitted noises
iN,   iB,     iK        = Noises index

This function must receive an Allan deviation in the standard format,
i.e., just one inflexion point, check mark shape.

ref: https://www.mathworks.com/help/fusion/ug/inertial-sensor-noise-analysis-using-allan-variance.html
"""
function allan_fit_simple(taus, adev)

    x = log10.(taus)
    y = log10.(adev)
    dy = diff(y) ./ diff(x)

    # ARW
    slope = -1/2 
    tau = 1                         # 1 s
    iN = argmin(abs.(dy.-slope))    # finding the index of minimum local error slope
    b = y[iN] - slope*x[iN]         # fitting a line with slope -1/2 
    logN = slope*log(tau) + b       # finding the value of the line at tau=1
    arw = 10^logN               

    # Bias drift
    slope = 0
    iB = argmin(abs.(dy.-slope))
    b = y[iB] - slope*x[iB]
    scfB = sqrt(2*log(2)/pi)  # 0.6643
    logB = b - log10(scfB)
    bias = 10^logB
    # If it does not work, try the minimum value instead.    
    # (val, iB) = findmin(adev)
    # bias = val/scfB

    # RRW
    slope = 1/2 
    tau = 3   # 3s or at 3h
    iK = argmin(abs.(dy.-slope))
    b = y[iK] - slope*x[iK]
    logK = slope*log10(tau) + b
    rrw = 10^logK

    return (arw=arw, bias=bias, rrw=rrw, ind_arw=iN, ind_bias=iB, ind_rrw=iK)
end
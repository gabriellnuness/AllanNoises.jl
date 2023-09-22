
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

    Q    = (avar_fit[1]/3)^(0.5)    # C = 3*Q^2; Q^2 = (q^2)/12, where q = quantum step size
    arw  = avar_fit[2]^(0.5)                        # C = N^2
    bias = (avar_fit[3]^(0.5)) / sqrt(2*log(2)/pi)  # C = (0.6643*B)^2
    rrw  = (3*avar_fit[4])^(0.5)                    # C = K^2 /3
    rr   = (2*avar_fit[5])^(0.5)                    # C = R^2 / 2

    avar_fit = sqrt.(avar_fit'*T')
    avar_fit = avar_fit[1,:]

    return (avar_fit, Q, arw, bias, rrw, rr)
end

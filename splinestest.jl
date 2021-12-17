using Plots
using Distributions, Random
using DSP
using LinearAlgebra

#Plot Styling
gr(size = (1000, 300), legend = false) 


#Set up true time series
t0 = 0
tn = 16

function target(x)
    return(sin(0.1*x^2) +sin(x) +  cos(0.08*x)^3 + 0.3*sin(x)^2)
end

#Set up x-axis
x = t0:0.1:tn

#Add noise to create target function
yTrue = target.(x)
yNoise = target.(x) + rand(Normal(0,0.5), length(x))

t = collect(x)
data = yNoise

#########Solution Starts Here


function pm(a)
    if a > 0
        return a
    end
    if a <= 0
        return 0
    end
end


function DesignMatrixCubic(t)
    B0 = ones(length(t))
    B1 = t
    B2 = t.^2
    B3 = t.^3
    return(hcat(B0,B1,B2,B3))
end


function TruncatedPowerBasisLinear(t,y)
    knots = [4,8,12]
    B0 = ones(length(t))
    B1 = t
    B2 = pm.(t .- knots[1])
    B3 = pm.(t .- knots[2])
    B4 = pm.(t .- knots[3])

    X = hcat(B0,B1,B2,B3,B4)

    coe = (X' * X *1)\(X'y*1)

    return(coe[1] .+ coe[2]*t + coe[3]*pm.(t .- knots[1]) + coe[4] *pm.(t .- knots[2])  + coe[5] *pm.(t .- knots[3]))

end

function TruncatedPowerBasisCubic(t,y)
    knots = [4,8,12]
    B0 = ones(length(t))
    B1 = t
    B2 = t.^2
    B3 = t.^3
    B4 = pm.(t .- knots[1]).^3
    B5 = pm.(t .- knots[2]).^3
    B6 = pm.(t .- knots[3]).^3

    X = hcat(B0,B1,B2,B3,B4,B5,B6)

    coe = (X' * X *1)\(X'y*1)

    return(coe[1] .+ coe[2]*t + coe[3]*t.^2 + coe[4]*t.^3 +
           coe[5]*pm.(t .- knots[1]).^3 + coe[6] *pm.(t .- knots[2]).^3  + coe[7] *pm.(t .- knots[3]).^3)

end



soln  = TruncatedPowerBasisLinear(t,yNoise)

soln  = TruncatedPowerBasisCubic(t,yNoise)


plot(t,yNoise, seriestype = :scatter, label = "data", markershape = :x, grid=false, ticks = false, showaxis = false);
hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);
plot!(t,yTrue, label = "true function", linestyle = :dash, linecolor = :black, linealpha = 1);
plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1)
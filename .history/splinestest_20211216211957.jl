using Plots
using Distributions, Random
using DSP

t0 = 0
tn = 16

function target(x)
    return(sin(0.1*x^2) +sin(x) +  cos(0.08*x)^3 + 0.3*sin(x)^2)
end

function cubic(t, x0, x1, x2, x3)
    return(x0 .+ x1.*t +x2.*t.^2 + x3.*t.^3)
end

function SubPlus(a)
    if a > 0
        return(a)
    end
    if a <= 0
        return(0)
    end
end


function MakeDesignMatrix(t)
    B1 = ones(length(t))
    B2 = t
    B3 = (SubPlus.(t .- 8).^3 .- SubPlus.(t .- t[end]).^3) ./ (t[end] .- 8)
    return(hcat(B1,B2,B3))
end


function FitSpline(X,y)
    coeff = inv(X'*X)*X'*y
    return(coeff)
end

function CubicFit(t,y)
    X = hcat(ones(length(t)), t,t.^2, t.^3);
    coeff = inv(X'*X)*X'*y;
    return(coeff) 
end

x = t0:0.1:tn

knots = t0:4:tn;

function pwCubic(t,y,knots)
    basisVector = []
    for i in knots
        y = cubicFit(t,y)
        hcat(basisVector,y)
    end
end




yTrue = target.(x)
yNoise = target.(x) + rand(Normal(0,0.5), length(x))


t = collect(x)
data = yNoise

cubicFit(t,data)

coeff = FitSpline(q,yNoise)

yPred = coeff[1] .+ coeff[2]*t .+ coeff[3]*t.^2 .+ coeff[4]*t.^3

yPred = cubic(t,coeff[1],coeff[2],coeff[3],coeff[4])

plot(x,yNoise, seriestype = scatter);
plot!(x,yTrue);
plot!(x,yPred)


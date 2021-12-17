using Plots
using Distributions, Random
using DSP

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


function cubic(t, x0, x1, x2, x3)
    return(x0 .+ x1.*t +x2.*t.^2 + x3.*t.^3)
end

plot(t,yNoise, seriestype = :scatter, label = "data", markershape = :x, grid=false, ticks = false, showaxis = false);
hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);
plot!(t,yTrue, label = "true function", linestyle = :dash, linecolor = :black, linealpha = 1)
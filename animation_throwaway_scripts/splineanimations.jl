##This is for the README header image

##Really rough code written very quickly using multiple cursors a lot so I don't need to think about functions. 
###Setup
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
y = yNoise

t = collect(x)
data = yNoise

##Make Cubic Spline Animation

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

soln = coe[1] .+ coe[2]*t + coe[3]*t.^2 + coe[4]*t.^3 +
       coe[5]*pm.(t .- knots[1]).^3 + coe[6] *pm.(t .- knots[2]).^3  + coe[7] *pm.(t .- knots[3]).^3



plot(t,yNoise, seriestype = :scatter, label = "data", markershape = :x, grid=false, ticks = false, showaxis = false);
hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);
plot!(t,yTrue, label = "true function", linestyle = :dash, linecolor = :black, linealpha = 1);
plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1);



@gif for i in reverse(10 .^(range(-6,stop=0.1,length=100)))
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1].+ i) * B0 
    b = (coe[2].+ i) * B1 
    c = (coe[3].+ i) * B2  
    d = (coe[4].+ i) * B3  
    e = (coe[5].+ i) * B4  
    f = (coe[6].+ i) * B5 
    g = (coe[6].+ i) * B6 

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end




##BASIS SPLINES

function pm2(x, a,b)
    if x >= a && x < b
        return(1)
    else
        return(0)
    end
end

knots = [2.8, 5.6, 8.4, 11.2, 14]
x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
y = [1,4,21,18,17,22,15,14,10,9,11,8,13,13,17]
#Pad n times
knots = [0, 0, 0, 0, 2.8, 5.6, 8.4, 11.2, 14, 15, 15, 15, 15]
B00 = pm2.(x, knots[1],knots[2])
B10 = pm2.(x, knots[2],knots[3])
B20 = pm2.(x, knots[3],knots[4])
B30 = pm2.(x, knots[4],knots[5])
B40 = pm2.(x, knots[5],knots[6])
B50 = pm2.(x, knots[6],knots[7])
B60 = pm2.(x, knots[7],knots[8])
B70 = pm2.(x, knots[8],knots[9])
B80 = pm2.(x, knots[9],knots[10])
B90 = pm2.(x, knots[10],knots[11])


B01 = replace!(((x .- knots[1])./(knots[1+1] - knots[1])).*B00, NaN =>0) + replace!(B10.*((knots[1+2] .- x)/(knots[1+2] - knots[1+1])), NaN => 0)
B11 = replace!(((x .- knots[2])./(knots[2+1] - knots[2])).*B10, NaN =>0) + replace!(B20.*((knots[2+2] .- x)/(knots[2+2] - knots[2+1])), NaN => 0)
B21 = replace!(((x .- knots[3])./(knots[3+1] - knots[3])).*B20, NaN =>0) + replace!(B30.*((knots[3+2] .- x)/(knots[3+2] - knots[3+1])), NaN => 0)
B31 = replace!(((x .- knots[4])./(knots[4+1] - knots[4])).*B30, NaN =>0) + replace!(B40.*((knots[4+2] .- x)/(knots[4+2] - knots[4+1])), NaN => 0)
B41 = replace!(((x .- knots[5])./(knots[5+1] - knots[5])).*B40, NaN =>0) + replace!(B50.*((knots[5+2] .- x)/(knots[5+2] - knots[5+1])), NaN => 0)
B51 = replace!(((x .- knots[6])./(knots[6+1] - knots[6])).*B50, NaN =>0) + replace!(B60.*((knots[6+2] .- x)/(knots[6+2] - knots[6+1])), NaN => 0)
B61 = replace!(((x .- knots[7])./(knots[7+1] - knots[7])).*B60, NaN =>0) + replace!(B70.*((knots[7+2] .- x)/(knots[7+2] - knots[7+1])), NaN => 0)
B71 = replace!(((x .- knots[8])./(knots[8+1] - knots[8])).*B70, NaN =>0) + replace!(B80.*((knots[8+2] .- x)/(knots[8+2] - knots[8+1])), NaN => 0)
B81 = replace!(((x .- knots[9])./(knots[9+1] - knots[9])).*B80, NaN =>0) + replace!(B90.*((knots[9+2] .- x)/(knots[9+2] - knots[9+1])), NaN => 0)


B02 = replace!(((x .- knots[1])./(knots[1+2] - knots[1])).*B01, NaN => 0) + replace!(B11.*((knots[1+3] .- x)/(knots[1+3] - knots[1+1])),NaN => 0)
B12 = replace!(((x .- knots[2])./(knots[2+2] - knots[2])).*B11, NaN => 0) + replace!(B21.*((knots[2+3] .- x)/(knots[2+3] - knots[2+1])),NaN => 0)
B22 = replace!(((x .- knots[3])./(knots[3+2] - knots[3])).*B21, NaN => 0) + replace!(B31.*((knots[3+3] .- x)/(knots[3+3] - knots[3+1])),NaN => 0)
B32 = replace!(((x .- knots[4])./(knots[4+2] - knots[4])).*B31, NaN => 0) + replace!(B41.*((knots[4+3] .- x)/(knots[4+3] - knots[4+1])),NaN => 0)
B42 = replace!(((x .- knots[5])./(knots[5+2] - knots[5])).*B41, NaN => 0) + replace!(B51.*((knots[5+3] .- x)/(knots[5+3] - knots[5+1])),NaN => 0)
B52 = replace!(((x .- knots[6])./(knots[6+2] - knots[6])).*B51, NaN => 0) + replace!(B61.*((knots[6+3] .- x)/(knots[6+3] - knots[6+1])),NaN => 0)
B62 = replace!(((x .- knots[7])./(knots[7+2] - knots[7])).*B61, NaN => 0) + replace!(B71.*((knots[7+3] .- x)/(knots[7+3] - knots[7+1])),NaN => 0)
B72 = replace!(((x .- knots[8])./(knots[8+2] - knots[8])).*B71, NaN => 0) + replace!(B81.*((knots[8+3] .- x)/(knots[7+3] - knots[7+1])),NaN => 0)

B03 = replace!(((x .- knots[1])./(knots[1+3] - knots[1])).*B02, NaN => 0) + replace!(B12.*((knots[1+4] .- x)/(knots[1+4] - knots[1+1])),NaN => 0)
B13 = replace!(((x .- knots[2])./(knots[2+3] - knots[2])).*B12, NaN => 0) + replace!(B22.*((knots[2+4] .- x)/(knots[2+4] - knots[2+1])),NaN => 0)
B23 = replace!(((x .- knots[3])./(knots[3+3] - knots[3])).*B22, NaN => 0) + replace!(B32.*((knots[3+4] .- x)/(knots[3+4] - knots[3+1])),NaN => 0)
B33 = replace!(((x .- knots[4])./(knots[4+3] - knots[4])).*B32, NaN => 0) + replace!(B42.*((knots[4+4] .- x)/(knots[4+4] - knots[4+1])),NaN => 0)
B43 = replace!(((x .- knots[5])./(knots[5+3] - knots[5])).*B42, NaN => 0) + replace!(B52.*((knots[5+4] .- x)/(knots[5+4] - knots[5+1])),NaN => 0)
B53 = replace!(((x .- knots[6])./(knots[6+3] - knots[6])).*B52, NaN => 0) + replace!(B62.*((knots[6+4] .- x)/(knots[6+4] - knots[6+1])),NaN => 0)
B63 = replace!(((x .- knots[7])./(knots[7+3] - knots[7])).*B62, NaN => 0) + replace!(B72.*((knots[7+4] .- x)/(knots[7+4] - knots[6+1])),NaN => 0)

B63[end] = 1



X = hcat(B03,B13,B23,B33,B43,B53,B63)

coe = (X' * X *1)\(X'y*1)

soln = coe[1]*B03 + coe[2]*B13 + coe[3]*B23 + coe[4]*B33 + coe[5]*B43 + coe[6]*B53 + coe[7]*B63

q = reverse(10 .^range(-2.5,stop=0.1,length=100))
@gif for i in q
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1] - 2*i ) * B03 
    b = (coe[2] / (1 + q[1]*abs(coe[2])) ) * B13 
    c = (coe[3] / (1 + q[1]*abs(coe[3]))) * B23  
    d = (coe[4] / (1 + q[1]*abs(coe[4]))) * B33  
    e = (coe[5] / (1 + q[1]*abs(coe[5]))) * B43  
    f = (coe[6] / (1 + q[1]*abs(coe[6]))) * B53 
    g = (coe[7] / (1 + q[1]*abs(coe[7]))) * B63 

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end


@gif for i in q
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1] ) * B03 
    b = (coe[2] / (1 + i*abs(coe[2]))) * B13 
    c = (coe[3] / (1 + q[1]*abs(coe[3]))) * B23  
    d = (coe[4] / (1 + q[1]*abs(coe[4]))) * B33  
    e = (coe[5] / (1 + q[1]*abs(coe[5]))) * B43  
    f = (coe[6] / (1 + q[1]*abs(coe[6]))) * B53 
    g = (coe[7] / (1 + q[1]*abs(coe[7]))) * B63  

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end

@gif for i in q
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1] ) * B03 
    b = (coe[2] ) * B13 
    c = (coe[3] / (1 + i*abs(coe[3]))) * B23  
    d = (coe[4] / (1 + q[1]*abs(coe[4]))) * B33  
    e = (coe[5] / (1 + q[1]*abs(coe[5]))) * B43  
    f = (coe[6] / (1 + q[1]*abs(coe[6]))) * B53 
    g = (coe[7] / (1 + q[1]*abs(coe[7]))) * B63 

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end
 
@gif for i in q
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1] ) * B03 
    b = (coe[2] ) * B13 
    c = (coe[3] ) * B23  
    d = (coe[4] / (1 + i*abs(coe[4]))) * B33  
    e = (coe[5] / (1 + q[1]*abs(coe[5]))) * B43  
    f = (coe[6] / (1 + q[1]*abs(coe[6]))) * B53 
    g = (coe[7] / (1 + q[1]*abs(coe[7]))) * B63 

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end
@gif for i in q
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1] ) * B03 
    b = (coe[2] ) * B13 
    c = (coe[3] ) * B23  
    d = (coe[4] ) * B33  
    e = (coe[5] / (1 + i*abs(coe[5]))) * B43  
    f = (coe[6] / (1 + q[1]*abs(coe[6]))) * B53 
    g = (coe[7] / (1 + q[1]*abs(coe[7]))) * B63 

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end

@gif for i in q
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1] ) * B03 
    b = (coe[2] ) * B13 
    c = (coe[3] ) * B23  
    d = (coe[4] ) * B33  
    e = (coe[5] ) * B43  
    f = (coe[6] / (1 + i*abs(coe[6]))) * B53 
    g = (coe[7] / (1 + q[1]*abs(coe[7]))) * B63 

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end

@gif for i in q
    plot(t,yNoise, seriestype = :scatter, label = "data",markershape = :x, markercolor = :black, grid=false, ticks = false, showaxis = false);
    #hline!([0],linestyle = :solid, linecolor = :grey, linealpha = 0.5);

    a = (coe[1] ) * B03 
    b = (coe[2] ) * B13 
    c = (coe[3] ) * B23  
    d = (coe[4] ) * B33  
    e = (coe[5] ) * B43  
    f = (coe[6] ) * B53 
    g = (coe[7] / (1 + i*abs(coe[7]))) * B63 

    plot!(t,a,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,b,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,c,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,d,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,e,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,f,linealpha = 0.5, palette = :Dark2_5);
    plot!(t,g,linealpha = 0.5, palette = :Dark2_5);

    soln = a + b + c +d + e + f + g
    plot!(t,soln, label = "estimation", linestyle = :solid, linecolor = :red, linealpha = 1, linewidth = 1.5, ylims=(floor(minimum(yNoise)),ceil(maximum(yNoise))))  

end
 
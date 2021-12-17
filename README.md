# Cpp Splines

![header](https://raw.githubusercontent.com/JasonPekos/CppSplines/main/cubicanim.gif)

## Scope
This repo contains object oriented spline-fitting code written in c++. Given some time series data, and command line arguments for:

- Number of knot points
- Class of basis functions
    - Truncated power basis
    - B-spline basis
    - ?
- ?

Following e.g. the SKlean style of statistics model classes, we have a model *(a)* instantiated with hyperparameters, *(b)* with .fit and  *(c)* with .predict methods. 

Return the corresponding fit. Plotting code is also provided. 
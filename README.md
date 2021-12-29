# Cpp Splines

![header](https://raw.githubusercontent.com/JasonPekos/CppSplines/main/cubicanim.gif)

## Scope
This repo contains object oriented spline-fitting code written in c++. Given some time series data, and command line arguments for:


- Class of basis functions
    - Truncated power basis
    - B-spline basis
- Degree of basis function
- Number of knot points


Following e.g. the SKlean style of statistics model classes, we have a model instantiated with hyperparameters, using methods to fit to data and to predict outcomes. 

Plotting code is also provided. 


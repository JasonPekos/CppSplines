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

**Overview:**

Although linear models are often convenient (or even necessary), in realistic applications, whatever function $f(x)$ we are estimating is highly likely to be non-linear. This package contains some basis expansion methods for fitting non-linear models. 

The common theme is that, instead of regressing on our vector of inputs directly, we instead augment or replace each input with a corresponding nonlinear transformation at that point, such that the transformations sum linearly to model target output. 

That is, instead of modelling the relationship as:

$$f(x) = \beta_0 + \beta_1x$$

We construct $M$ nonlinear transformations $b_m(x)$ for each $x$ and $m \leq M$, and model:

$$f(x) = \sum_0^M \beta_m (b_m(x))$$

Where the utility of this approach comes from the fact that, after applying the non-linear transformation, the functions enter the model linearly. This means that the vast array of linear model literature can — in many cases — be adapted to apply to basis expansion methods, a benefit which doesn't appear when using e.g. gradient boosted trees or deep neural networks. 

Additionally, unlike some other machine learning methods, these remain (weakly) convex optimization problems, and so we are assured that, if a solution is recovered, it is a global maxima. 

**Regression Splines**

The first method available here is a regression spline. We split the data along the time axis into $k$ interior knots, and fit a piecewise polynomial between each knot. Additionally, we enforce smoothness constraints such that each successive polynomial has continuous value, first, and second derivatives at each knot point where it meets. We use two different basis construction methods here --- BSplines and Power Basis splines. The BSPlines represent a natural spline basis, in contrast to the Power Basis, which enforces no such boundary constraint. 

*Power Basis*

The power basis representation of a spline is found by augmenting the basis for a polynomial regression problem of that same order with additional basis functions at each not to enforce our continuity constraints.

For some degree $n$ basis with $k$ knots, our basis representation is given by:


$$\begin{array}{lll}
b_{1}(x)=1, & b_{2}(x)=x, & b_{3}(x)=x^{2}, & b_{n}(x) = x^n \\
\end{array}$$

$$\begin{array}{lll}
b_{n+1}(x) = \left(x-k_{1} \right)_{+}^{n} , & b_{n+2}(x) = \left(x-k_{2} \right)_{+}^{n}  \\
\end{array}$$

$$\begin{array}{lll}
b_{n+k}(x) = \left(x-k_{k} \right)_{+}^{n}  \\
\end{array}$$

Where $(x)_+$ is zero if $x$ is negative, and $x$ otherwise. 

For data in 'input.csv', we can return a power basis regression spline with:

```{bash}
./splines PowerBasis power knots
```

for example, 

```{bash}
./splines PowerBasis 3 2
```

Returns a regression spline with degree $3$ with $2$ knots. The output is given as time series data in 'output.csv'. Example output for the above code is given by:

| t | y |
| --- | ---|
|  1    | -0.947246   |
|  1.1    |  0.320104  |
|  1.2    |  1.53814  |
|  1.3    |  2.70769  |
|   ...   |  ...  |

Which can be automatically plotted with the attached plotting.py file, returning:

![plotone](https://raw.githubusercontent.com/JasonPekos/CppSplines/main/images/PowerBasis32.png)

(using the default input.csv provided as data).

Notes:

- Asking for more knots than datapoints probably won't, for nearly all data structure. 

- Longer computations should use B-Splines instead. 

*B-Splines*

Although computationally simple, the power basis provided above has numerous undesirable numerical properties --- specifically, they require the computation of large numbers, which can lead to rounding problems. To solve this issue, we turn towards an alternative basis construction which can span the same set of functions as a clamped power basis. 

Basis spline are wonderful in that they provide local support, and are much more well behaved numerically. 

Unlike the power basis, higher order basis functions have no simple closed form, and must be represented as the iterative convolution of lower order basis functions, where the $0$th order is an indicator function on successive knots.

For data in 'input.csv', we can return a power B-spline regression spline with:

```{bash}
./splines BSpline power knots
```

for example, 

```{bash}
./splines BSpline 3 6
```

Returns a regression spline with degree $3$ with $6$ interior knots. The output is given as time series data in 'output.csv'. Example output for the above code is given by:

| t | y |
| --- | ---|
|  1    | 0.993459   |
|  1.1    |  -1.38922  |
|  1.2    |  -2.93915  |
|  1.3    | -3.74152  |
|   ...   |  ...  |

Which can be automatically plotted with the attached plotting.py file, returning:


![plotone](https://raw.githubusercontent.com/JasonPekos/CppSplines/main/images/BSpline36.png)



**Smoothing Splines**














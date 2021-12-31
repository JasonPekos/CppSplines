#include <vector>
#include <cmath>

#pragma once

double pm(double x)
{
    /**
     * @brief This returns 'x' iff x > 0. 
     * 
     * Used often in power spline basis construction. 
     */
    if (x > 0)
    {
        return (x);
    }
    else
    {
        return (0);
    }
}

double CoxDeBoor(double x, uint64_t index, std::vector<double> knots, uint64_t power)
{
    /**
     * @brief Compute the basis spline at 'x' at knot 'index' recursively using De Boor's recursion formula.
     * 
     * @param x B_nm(x); basis function at this point
     * @param index Knot index
     * @param knots Vector of knot positions
     * @param power Power of B-Spline
     * 
     * @return the basis with the above parameters evaluated at x.
     */

    //In the zero degree case, we simply have a step function.
    double out = 0;
    double a1 = 0;
    double a2 = 0;

    //Quick check to avoid division by zero.
    if (power == 0)
    {
        if (x >= knots[index])
        {
            /*
            These constitute the base case in a recursive function.

            Note: this is simply a step function, the 0th order BSpline basis. 
            */
            if (x < knots[index + 1])
            {
                out = 1;
                return (out);
            }
            else
            {
                out = 0;
                return (out);
            }
        }
        else
        {
            out = 0;
            return (out);
        }
    }

    //Define 0/0 := 0 for this problem.
    if (knots[power + index] == knots[index])
    {
        a1 = 0;
    }
    else //If no division by zero, apply Cox De Boor formula
    {
        a1 = (x - knots[index]) / (knots[power + index] - knots[index]);
    }

    //Avoiding another division by zero error.
    if (knots[power + index + 1] == knots[index + 1])
    {
        a2 = 0;
    }
    else
    {
        a2 = (knots[index + power + 1] - x) / (knots[power + index + 1] - knots[index + 1]);
    }

    //Apply this recursively, wrapping up from the step function base case up to the 'power' degree basis spline.
    out = a1 * CoxDeBoor(x, index, knots, power - 1) + a2 * CoxDeBoor(x, index + 1, knots, power - 1);

    return (out);
}

std::vector<double> linspace(double a, double b, double by)
{
    /**
     * @brief linearly spaced points between a,b with spacing 'by'.
     * 
     * @param a start point
     * @param b end point
     * @param by spacing between points. 
     * 
     */
    std::vector<double> out = {};
    double current = a;

    while (current < b)
    {
        out.push_back(current);
        current += by;
    }

    return (out);
}
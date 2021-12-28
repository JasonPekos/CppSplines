#include <vector>
#include <cmath>

#pragma once

double pm(double x){
    /**
     * @brief This returns 'x' iff x > 0. 
     * 
     * Used often in power spline basis construction. 
     */
    if (x > 0)
    {
        return(x);
    }
    else
    {
        return(0);
    }
}

double CoxDeBoor(double x, uint64_t index, std::vector<double> knots, uint64_t power){
    /**
     * @brief Computer the basis spline at x "at" knot 'index' recursively using De Boor's recursion formula (DeBoor 2001)
     * 
     * @param x B_nm(x); basis function at this point
     * @param index Knot index
     * @param knots Vector of knot positions
     * @param power Power of B-Spline
     */

    //In the zero degree case, we simply have a step function. 
    double out = 0;
    double a1  = 0;
    double a2  = 0;
    if (power == 0)
    { 
        if (x >= knots[index])
        {
            if (x < knots[index+1])
            {
                out = 1;
                return(out);
            }
            else
            {
                return(0);
            }
        }
        else
        {
            return(0);
        }
    }
    if ((knots[power + index] - knots[index]) == 0)
    {
        a1 = 0; 
    }
    else
    {
        a1 = (x - knots[index])/(knots[power + index] - knots[index]);
    }
    if ((knots[power + index + 1] - knots[index + 1]) == 0)
    {
        a2 = 0;
    }
    else
    {
        a2 = (knots[index + power + 1] - x)/(knots[power + index + 1] - knots[index + 1]);
    }

    out = a1*CoxDeBoor(x,index, knots, power-1) + a2*CoxDeBoor(x,index + 1, knots, power -1);
    
    return(out);
}
#include <iostream>
#include <vector>
using namespace std;
#include "linAlgHelperFunctions.hpp" 
#pragma once

//Ben Klemens:

/*
X -> Beta (given data, fit model)
Beta -> X (given model, generate data)
(X, B) -> (Given both data and parameters, estimate their likelihood or probability)

Focusing mainly on first option here re: model functionality.   
*/


class Spline
{
private:

    uint64_t Knots = 3; //Data about model (number of knots)
    uint64_t Power = 3; //Data about model (power of knots)
    std::string Method = "PowerBasis"; //Data about model (Method) [PowerBasis or BSpline or PolynomialRegression] 
    
    
public:

std::vector<std::vector<double>> Coe{1,std::vector<double>(1,0)};

void fit(std::vector<double> t, std::vector<double> y){
        //estimate(data[0], data[1], Coe (by reference ))

        //If asked for PowerBasis
        if (Method == "PowerBasis")
        {
            std::cout << "Not Implemented Yet" << "\n";
        }

        //If asked for Basis Spline
        if (Method == "BSpline")
        {
           std::cout << "Not Implemented Yet" << "\n";
        }

        //If asked for polynomial regression. 
        if (Method == "PolynomialRegression")
        {
            if (Knots > 0)
            {
                std::cout << "Warning: Knots parameter only applies to Bspline and PowerBasis methods.";
            }
            std::vector<std::vector<double>> DesignMatrix = Design(t,Power);
            Coe = SolveSystem(MatMul(Transpose(DesignMatrix), DesignMatrix), MatVecMul(Transpose(DesignMatrix), y));
        }
    }

    double predict(double t){
        double val = 0;
        if (Method == "PolynomialRegression")
        {
            
            for (uint64_t i = 0; i < Coe.size(); i++)
            {
                val += Coe[i][0]*pow(t, i);
            }
        }
        return(val);
    }
    Spline(std::string method, uint64_t power, uint64_t knots){
        /* constructor --- add error checking 

        knots -> int64_t, positive value
        method -> string, matches "PowerBasis, BSpline"
        power -> int64_t, nth order splines, positive or zero, 


        */
        Knots = knots;
        Power = power;
        Method = method;
    }
};




#include <iostream>
#include <vector>
using namespace std;
#include "linAlgHelperFunctions.hpp" 
#include "mathHelperFunctions.hpp"
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
     //Holds the actual knot positions. 

public:

std::vector<std::vector<double>> Coe{1,std::vector<double>(1,0)};
std::vector<double> kTemp = {};

void fit(std::vector<double> t, std::vector<double> y){

        //If asked for PowerBasis
        if (Method == "PowerBasis")
        {
            kTemp = {};
            if (Knots > 0)
            {
                double range = t.back() - t.front();
                double kEvery  = range / Knots;
                double kCurrent = kEvery;
                while (kCurrent < t.back())
                {
                    kTemp.push_back(kCurrent);
                    kCurrent += kEvery;
                }
                
            }
            else
            {
                std::cout << "error: at least one knot required" << "\n";
                exit(-1);
            }

            std::vector<std::vector<double>> DesignMatrix = DesignPowerBasis(t, Power, kTemp);


            std::vector<std::vector<double>> XTX = MatMul(Transpose(DesignMatrix), DesignMatrix);

            Coe = SolveSystem(XTX, MatVecMul(Transpose(DesignMatrix), y));
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
        if (Method == "PowerBasis")
        {
            //Before the power basis functions. 

            val = 0;
            for (uint64_t i = 0; i < Power + 1; i++)
            {
                val += Coe[i][0]*pow(t, i);
            }
            
            for (uint64_t i = Power + 1; i < Coe.size(); i++)
            {
                val += Coe[i][0]*pow(pm(t - kTemp[i - (Power + 1 )]), Power);
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




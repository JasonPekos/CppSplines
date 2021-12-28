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

public:

//Holds the coefficients. Initialized at zero.
std::vector<std::vector<double>> Coe{1,std::vector<double>(1,0)};

//Holds the actual knot positions. Filled when fitting the model.
std::vector<double> kTemp = {};

void fit(std::vector<double> t, std::vector<double> y){
    /**
     * @brief Fit the model on data (t,y).
     * 
     */

        //If asked for PowerBasis, continue by fitting with the following method:
        if (Method == "PowerBasis")
        {
            //Empty out the knot positions vector, so we can push into the back safely. 
            kTemp = {}; 
            if (Knots > 0)
            {
                /*
                Place 'Knots' evenly spaced knots over the range of time values. 
                */
                double range = t.back() - t.front();
                double kEvery  = range / (Knots + 1);
                double kCurrent = kEvery;
                while (kCurrent < t.back())
                {
                    kTemp.push_back(kCurrent);
                    kCurrent += kEvery;
                }
                
            }
            else
            {
                /*
                If this is somehow being called without CL args (e.g. during code reuse), call an error if knots <= 0. 
                */
                std::cout << "error: at least one knot required" << "\n";
                exit(-1);
            }

            //Create the design matrix for this program using the LinAlg function. 
            std::vector<std::vector<double>> DesignMatrix = DesignPowerBasis(t, Power, kTemp);

            //Solve as a linear system. 
            std::vector<std::vector<double>> XTX = MatMul(Transpose(DesignMatrix), DesignMatrix);

            //Return Coefficients.
            Coe = SolveSystem(XTX, MatVecMul(Transpose(DesignMatrix), y));
        }



        //If asked for polynomial regression. 
        if (Method == "PolynomialRegression")
        {
            //Warning if you asked for knots on a method that doesn't support them.
            if (Knots > 0) 
            {
                std::cout << "Warning: Knots parameter only applies to Bspline and PowerBasis methods.";
            }
            //Create design matrix and solve as a linear system, e.g. X'X \ X'y
            std::vector<std::vector<double>> DesignMatrix = Design(t,Power);
            Coe = SolveSystem(MatMul(Transpose(DesignMatrix), DesignMatrix), MatVecMul(Transpose(DesignMatrix), y));
        }

        //If asked for Basis Spline
        if (Method == "BSpline")
        {
           std::cout << "Not Implemented Yet" << "\n";
        }
    }

    double predict(double t){
        /**
         * @brief Predict output given model + input
         * 
         * @param t Value we need to predict at, e.g. return Model(t)
         * 
         * -> sum(coe_n*basisFuction_n (t)) for all n. 
         * 
         */

        //Define value we sum into. 
        double val = 0;
        if (Method == "PolynomialRegression")
        {
            //Again, sum(coe_n*basisFuction_n (t)) for all n.
            for (uint64_t i = 0; i < Coe.size(); i++)
            {
                val += Coe[i][0]*pow(t, i);
            }
        }
        if (Method == "PowerBasis")
        {
            //Before the explicitly power basis functions. 
            val = 0;
            for (uint64_t i = 0; i < Power + 1; i++)
            {
                val += Coe[i][0]*pow(t, i);
            }
            //Power basis functions, e.g. coe * +_(t - knot)^max_power. 
            for (uint64_t i = Power + 1; i < Coe.size(); i++)
            {
                val += Coe[i][0]*pow(pm(t - kTemp[i - (Power + 1 )]), Power);
            }
        }
        return(val);
    }

    Spline(std::string method, uint64_t power, uint64_t knots){
        /**
         * @brief constructor
         * 
         * @param knots int64_t, positive value.
         * @param method string, matches "PowerBasis, BSpline", "PolynomialRegression".
         * @param power int64_t, nth order splines, positive or zero.
         */

        //Strings to check against in constructor. 
        std::string BSP = "BSpline";
        std::string PBS = "PowerBasis";
        std::string PR  = "PolynomialRegression";

        //Assign values. 
        Knots = knots;
        Power = power;
        Method = method;

        //Simple IO check in the case of explicit construction outside of CL arguments --- makes the class safer. 
        if (knots < 0)
        {
            std::cout << "Error: can't construct with a negative number of knots!" << "\n";
            throw std::invalid_argument( "received negative value" );
        }
        if (power < 1)
        {
            std::cout << "Error: can't construct with a power less than one!" << "\n";
            throw std::invalid_argument( "received illegitimate value" );
        }
        if (method != BSP && method != PBS && method != PR )
        {
            std::cout << "Error: submit a valid method!" << "\n";
            throw std::invalid_argument( "received illegitimate method" );
        }
    }
};




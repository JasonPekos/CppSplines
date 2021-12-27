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

    int64_t Knots = 3; //Data about model (number of knots)
    int64_t Power = 3; //Data about model (power of knots)
    std::string Method = "PowerBasis"; //Data about model (Method) [PowerBasis or BSpline] 

    std::vector<std::vector<double>> Coe(A.size(), std::vector<double>(1, 0));


    void fit(){

    }

    void predict(){

    }
    
public:
    Spline(std::string method, int64_t power, int64_t knots){
        /* constructor --- add error checking 

        knots -> int64_t, positive value
        method -> string, matches "PowerBasis, BSpline"
        power -> int64_t, nth order splines, positive or zero, 


        */
        Knots = knots;
        Power = power;
        Method = method;
    }
    ~Spline();
};




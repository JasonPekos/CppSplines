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
     * Uses the closed form solutions: Coefficients = (XT*X)^-1 * XT*y
     * 
     * Solved using Gaussian Elimination with Back Substitution instead of a matrix inverse, 
     * as numerical computation of the inverse is ill-conditioned. 
     * 
     */

        if (Knots > (t.size() - 1))
        //Quick warning if Knots > Datapoints. 
        {
            std::cout << "WARNING: initializing a spline with more knots than datapoints will likely lead to perfect multicolinearity. ";
            std::cout << "It is extremely unlikely that this model was fit correctly. \n";

        }

        //If asked for PowerBasis, continue by fitting with the following method:
        if (Method == "PowerBasis")
        {
            //Empty out the knot positions vector, so we can push into the back safely. 
            kTemp = {}; 
            if (Knots == 1)
            {
                //If we have one knot only, place it in the middle of the data. 
                kTemp.push_back((t.back() - t.front() )/ 2);
            }
            
            if (Knots > 1)
            {
                /*
                Place 'Knots' evenly spaced knots over the range of time values. 
                */
                double range = t.back() - t.front();
                double kEvery  = range / (Knots);
                double kCurrent = kEvery;
                while (kCurrent < t.back())
                {
                    kTemp.push_back(kCurrent);
                    kCurrent += kEvery;
                }
                
            }
            if (Knots < 1)
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
            //Empty out the knot positions vector, so we can push into the back safely. 
            kTemp = {}; 
            if (Knots == 1)
            {
                //If there is one single knot, place it in the middle of the interval.
                kTemp.push_back((t.back() - t.front() )/ 2);
            }
            
            if (Knots > 1)
            {
                /*
                Place 'Knots' evenly spaced knots over the range of time values. 
                */
                double range = t.back() - t.front();
                double kEvery  = range / (Knots);
                double kCurrent = kEvery;
                while (kCurrent < t.back())
                {
                    kTemp.push_back(kCurrent);
                    kCurrent += kEvery;
                }
                
            }
            if (Knots < 1)
            {
                /*
                If this is somehow being called without CL args (e.g. during code reuse), call an error if knots <= 0. 
                */
                std::cout << "error: at least one knot required" << "\n";
                exit(-1);
            }

            //Push in padding knots for basis splines. 
            uint64_t tempcounter = 0;
            while (tempcounter <= Power + 1) //Add 'power + 1' padded knots on each side. 
            {
                kTemp.insert(kTemp.begin(), 0);
                kTemp.push_back(t.back());
                tempcounter +=1;
            }

            //Create the design matrix for this program using the LinAlg function. 
            std::vector<std::vector<double>> DesignMatrix = DesignBSplineBasis(t, Power, kTemp);

            //Solve as a linear system. 
            std::vector<std::vector<double>> XTX = MatMul(Transpose(DesignMatrix), DesignMatrix);

            //Return Coefficients.
            Coe = SolveSystem(XTX, MatVecMul(Transpose(DesignMatrix), y));
        }
        if (MatNoNAN(Coe) == 0)
        {
            /*
            Some combinations of data / basis representation / knot amount can lead to perfectly multicolinear columns,
            and therefore NaNs after we try to solve the system. The resulting system cannot be plotted, so we 
            throw a warning here, and let the user know model fitting has failed. 

            We still save the coefficient matrix because the position of the NaNs can be informative. 
            */
            std::cout << "\n \n WARNING: The combination of knots / degrees / input data submitted has resulted in numerical instability. "; 
            std::cout << "Model fitting has failed. Recommendation is to reduce knots / power, especially for small datasets. \n \n";
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
        if (Method == "BSpline")
        {
            val = 0;
            //Again, sum(coe_n*basisFuction_n (t)) for all n.
            for (uint64_t i = 0; i < Coe.size(); i++)
            {
                val += Coe[i][0]*CoxDeBoor(t,i + 1, kTemp, Power);
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
        if (method == PBS)
        {
            if (power > 4)
            {
                std::cout << "WARNING: Asking for a power basis with power > 4." << "\n";
                std::cout << "Powers this large result in numerically unstable regression matrices" << "\n";
            }
        }        
    }
};

class GAM
{
private:
    uint64_t Power  = 3;
    std::vector<double> kTemp = {};

public:
    double Lambda = 0;
    std::vector<std::vector<double>> Coe{1,std::vector<double>(1,0)};
    double CrossValScore = 100;

void fit(std::vector<double> t, std::vector<double> y){
    /**
     * @brief Fit the model on data (t,y).
     * 
     * Uses the closed form solutions: Coefficients = (XT*X + Lambda * Omega)^-1 * XT*y
     * 
     * Solved using Gaussian Elimination with Back Substitution instead of a matrix inverse, 
     * as numerical computation of the inverse is ill-conditioned. 
     * 
     */
   
        //Empty out the knot positions vector, so we can push into the back safely. 
        kTemp = {}; 

        //Add a knot at each datapoint except for the boundary. 
        for (uint64_t i = 1; i < t.size(); i++)
        {
            kTemp.push_back(t[i]);
        }

        //Push in padding knots for basis splines. 
        uint64_t tempcounter = 0;
        while (tempcounter <= Power + 1) //Add 'power + 1' padded knots on each side. 
        {
            kTemp.insert(kTemp.begin(), 0);
            kTemp.push_back(t.back());
            tempcounter +=1;
        }

        //Create the design matrix for this program using the LinAlg function. 
        std::vector<std::vector<double>> DesignMatrix = DesignBSplineBasis(t, Power, kTemp);

        //Multiply XTX
        std::vector<std::vector<double>> XTX = MatMul(Transpose(DesignMatrix), DesignMatrix);

        //Add penalty
        std::vector<std::vector<double>> P = AddWigglyPenalty(Lambda, XTX);

        //Return Coefficients.
        Coe = SolveSystem(P, MatVecMul(Transpose(DesignMatrix), y));


    
        if (MatNoNAN(Coe) == 0)
        {
            /*
            Some combinations of data / basis representation / knot amount can lead to perfectly multicolinear columns,
            and therefore NaNs after we try to solve the system. The resulting system cannot be plotted, so we 
            throw a warning here, and let the user know model fitting has failed. 

            We still save the coefficient matrix because the position of the NaNs can be informative. 
            */
            std::cout << "\n \n WARNING: The combination of knots / degrees / input data submitted has resulted in numerical instability. "; 
            std::cout << "Model fitting has failed. Recommendation is to reduce knots / power, especially for small datasets. \n \n";
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

        //Again, sum(coe_n*basisFuction_n (t)) for all n.
        for (uint64_t i = 0; i < Coe.size(); i++)
        {
            val += Coe[i][0]*CoxDeBoor(t,i + 1, kTemp, Power);
        }
        
        return(val);
    }

    void GCV(std::vector<double> t, std::vector<double> y){
        fit(t,y);

        double GCV = 0;
        std::vector<double> yDiff = {};
        std::vector<double> yHatDiff = {};

        //Numerator.
        for (uint64_t i = 0; i < y.size(); i++)
        {
            GCV += y.size()*pow((y[i] - predict(t[i])),2);
        }

      
        GCV = GCV / pow((y.size() - TraceHatMatrix(DesignBSplineBasis(t, Power, kTemp), Lambda)),2);

        CrossValScore = GCV;
    }

    void seekLambda(std::vector<double> t, std::vector<double> y){
        std::vector<double> tempLambda = {};

        std::vector<double> testpoints = linspace(0.1,10,0.1);

        for (double i = 0; i < testpoints.size(); i++)
        {
            Lambda = testpoints[i];
            GCV(t,y);
            tempLambda.push_back(CrossValScore);
            //std::cout << testpoints[i] << "| GCV: " << CrossValScore << "\n"; Uncomment to watch the CV score go down (good) and then up again.
        }

        std::vector<double>::iterator iter = std::min_element(std::begin(tempLambda), std::end(tempLambda));
        int64_t index = std::distance(std::begin(tempLambda), iter);
        Lambda = testpoints[index];
    }



    GAM(double lambda){
        /**
         * @brief Constructor for a GAM.
         * 
         */
        Lambda = lambda;
    }
};






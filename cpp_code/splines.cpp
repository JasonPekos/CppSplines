/**
 * @file splines.cpp
 * @author Jason Pekos (pekosj at mcmaster (dot) ca)
 * @brief Contains code for fitting numerous regression spline / GAM style models, with extensions into e.g. GCV for automatic parameter inference.
 * @version 0.01
 * @date 2021-12-31
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "model.hpp"                 // Contains sklearn style classes for Regression and Smoothing splines.
#include "linAlgHelperFunctions.hpp" // Helper functions for LinAlg related compuations, including stats specific computations (e.g. adding a penalty in a regression problem).
#include "iocheck.hpp"               // Input output checking
#include <vector>                    // Class used for linalg computations, dynamic arrays.
#include <string>                    // Holding IO strings.
#include <iostream>                  // command line input, write errors to command line.
#include <fstream>                   // Read, write to csv.
#include <algorithm>                 // Find First Not Of
#include <cmath>                     // Various math functions, e.g. abs.

int main(int argc, char const *argv[])
{

    // Quick IO check, moved into a header file.
    if (InputCheck(argc, argv) == 0)
    {
        std::cout << "Failed Input Output Check"
                  << "\n";
        return (-1);
    }

    // Turn "method" command line argument into string.
    std::string method(argv[1]);

    // Strings to check against.
    std::string PNR = "PolynomialRegression";
    std::string BSP = "BSpline";
    std::string PSP = "PowerBasis";
    std::string SMS = "Smooth";

    // Set defaults for variables to avoid memory errors.
    uint64_t knots = 0;
    uint64_t power = 3;
    double lambda = 2;

    // This is a flag for if 'auto' is called as the 'lambda' argument, which requires some special handling.
    bool autoFlag = 0;

    // Add parameters from command line arguments based on method.
    if (method == SMS)
    {
        // SMS = Smooth Spline / GAM style model.
        std::string Auto = "auto";
        std::string Arg2(argv[2]);

        // Turn on 'auto' flag in case where it was called.
        if (Arg2 == Auto)
        {
            autoFlag = 1;
        }
        else
        {
            // Store 'lambda' as a double (this is the wiggliness parameter).
            lambda = stod(Arg2);
        }
    }

    // If we aren't using a smoothing spline, we have a different command line argument.
    if (method != SMS)
    {
        // The second argument is always the power of the spline. We fix the GAM spline to be third order, which is relatively standard practice.
        // See: Wood (2006), Hastie and Tibshirani ESL chapter 5.
        std::string Arg2(argv[2]);
        char* end;
        power = strtoull(Arg2.c_str(),&end,10);
    }

    // If we aren't using a smoothing spline or Polynomial Regression, we need to manually specify knot number.
    if (method != PNR && method != SMS)
    {
        // Smoothing splines place a knot at every data point, and Polynomial Regression has no knots.
        std::string Arg3(argv[3]);
        char* end;
        knots = strtoull(Arg3.c_str(),&end,10);
    }

    /*
    Concluding command line input, move on to reading input.csv.
    */

    // Empty vectors to push .CSV input back into
    std::vector<double> t = {};
    std::vector<double> y = {};

    // temp variables to store IOStream values.
    std::string xVal;
    std::string yVal;

    /*
    Read in input.csv file
    */

    // Inputstream to the .csv we care about.
    std::ifstream data("input.csv");

    // Check to make sure file opened correctly.
    if (!data.is_open())
    {
        std::cout << "Error opening file!";
        return (-1);
    }

    // Skip header.
    std::getline(data, xVal, ',');
    std::getline(data, yVal, '\n');

    // This checks to make sure that the first line is a header, not a number.
    if (xVal.find_first_not_of("-1234567890.") == std::string::npos)
    {
        std::cout << "warning! check to make sure .csv has column names. "
                  << "\n";
    }
    if (yVal.find_first_not_of("-1234567890.") == std::string::npos)
    {
        std::cout << "warning! check to make sure .csv has column names. "
                  << "\n";
    }

    // Read in data until we hit the end of the file
    while (data.peek() != EOF)
    {

        // Push data into xVal, yVal with ',' and newline delimiters.
        std::getline(data, xVal, ',');
        std::getline(data, yVal, '\n');

        // Check for legitimate input. Must be float. Can be negative and spline is still well defined.
        if (xVal.find_first_not_of("-1234567890.") != std::string::npos)
        {
            std::cout << "Error with .csv input file on line containing:" << xVal << "\n";
            return (-1);
        }

        // Similar check for legitimate input.
        if (yVal.find_first_not_of("-1234567890.") != std::string::npos)
        {
            std::cout <<  "Error with .csv input file on line containing:" << yVal << "\n";
            return (-1);
        }

        // Push out value for this loop into the data vector.
        t.push_back(std::stod(xVal));
        y.push_back(std::stod(yVal));
    }

    /*
    Final check to make sure we have time series data, and no corrupt cells.
    */
    if (t.size() != y.size())
    {
        std::cout << "Error --- input data of different lengths";
        return (-1);
    }

    /*
    All IO checks complete, moving on to fitting the model and outputting into a .csv file.
    */

    // Temp vectors to hold predictions for us to push into the output.csv file.
    std::vector<double> xTemp = linspace(t[0], t.back(), fabs((t.back() - t[0]) / 10000)); //Spaced like this to avoid weird time series input over tiny ranges.

    std::vector<double> modelTemp = {}; //Model output.

    /*
    If method is NOT a smoothing spline / GAM, we use the same regression spline class, and handle e.g. the basis inside the class call.

    This is because, unlike with GAMs, this class can be written so that every method applies to any of Polynomial Regression / BSplines / Power Basis Splines.

    GAMS have additional methods, and additional requirements around fitting models when incorporating the penalization term.
    */
    if (method != SMS)
    {
        /*
        Call Spline (again, not GAM) class, instantiated with (method / power / knots)

        method = first command line argument.
        power  = degree of the basis.
        knots  = number of knots.

        */
        Spline model(method, power, knots);

        // fit the model to time series data t,y with a method call.
        model.fit(t, y);

        // push the predictions over the smooth, finer-grained interval into the output storage vector.
        for (uint64_t i = 0; i < xTemp.size(); i++)
        {
            modelTemp.push_back(model.predict(xTemp[i]));
        }
    }
    else
    {
        /*
        The GAM model is always called with 'Lambda'.

        Power: Always third degree splines, as is standard. 
        Knots: GAMS place knots at every predictor and avoid overfitting via penalization, so this is not needed.
        Lambda: Penalization term. 
        */
        GAM model(lambda);

        if (autoFlag == 1)
        {
            // If we asked for 'auto' as a command line argument, we update the Lambda parameter for our instance of the GAM model with this method.
            model.seekLambda(t, y);

            // This method can take a while, so we let the user know when we've recovered lambda.
            std::cout << "Recovered Lambda: " << model.printLambda() << "\n";
        }

        // Fit the model, push to output, same as with the regression spline.
        model.fit(t, y);
        for (uint64_t i = 0; i < xTemp.size(); i++)
        {
            modelTemp.push_back(model.predict(xTemp[i]));
        }
    }

    /* 
    FILE OUTPUT
    */
    std::ofstream output("output.csv"); // IO check.
    if (!output.is_open())
    {
        std::cout << "Error opening output file: input.csv probably missing from cpp_code directory.";
        return (-1);
    }

    // Output CSV header names.
    output << "t" << ',';
    output << "y";
    output << '\n';

    // Populate .CSV with model output over correct range.
    for (uint64_t i = 0; i < modelTemp.size(); i++)
    {
        output << xTemp[i] << ",";
        output << modelTemp[i] << '\n';
    }

    output.close();

    // Message to annouce completion of program.
    std::cout << "output.csv updated \n"; //

    return 0;
}

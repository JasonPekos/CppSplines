#include "model.hpp" //Leave out until LA is working
#include "linAlgHelperFunctions.hpp" // when I move everything over.
#include "iocheck.hpp" //iochecking
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>



int main(int argc, char const *argv[])
{

    //Quick IO check, moved into a header file. 
    if(InputCheck(argc, argv) == 0){
        std::cout << "Failed Input Output Check" << "\n";
        return(-1);
    }

    //Method CL argument into string. 
    std::string method(argv[1]);

    //In the case of no knot points being provided, simply use the un-augmented polynomial regression function. 
    if (InputCheck(argc, argv) == 2)
    {
        method = "PolynomialRegression"; //If we don't have any knot points, just use PNR as the default. 
    }

    //Bring over knots and power CL arguments.
    std::string Arg2(argv[2]);

    std::string PNR  = "PolynomialRegression";
    std::string BSP  = "BSpline";
    std::string PSP  = "PowerBasis";
    

    //Case to fixed width values. 
    uint64_t power = (uint64_t)stoi(Arg2);
    uint64_t knots = 0;


    if (method != PNR)
    {
        std::string Arg3(argv[3]);
        knots = (uint64_t)stoi(Arg3);
    }
    

    //Empty vectors to push .CSV input back into
    std::vector<double> t = {}; 
    std::vector<double> y = {};

    //temp variables to store IOStream values.  
    std::string xVal; 
    std::string yVal;

    /*
    Read in input.csv file
    */

    //Inputstream to the .csv we care about. 
    ifstream data("input.csv");

    //Check to make sure file opened correctly. 
    if (!data.is_open()){
        cout << "Error opening file!";
        return(-1);
    }

    //Read in data until we hit the end of the file. 
    while (data.peek()!= EOF)
    {

        //Push data into xVal, yVal with ',' and newline delimiters.
        std::getline(data, xVal, ','); 
        std::getline(data, yVal, '\n');


        //Check for legitimate input. Must be float. Can be negative and spline still well defined. 
        if (xVal.find_first_not_of("-1234567890.")!= std::string::npos)
        {
            std::cout << "Error with .csv input file";
            return(-1);
        }

        if (yVal.find_first_not_of("-1234567890.")!= std::string::npos)
        {
            std::cout << "Error with .csv input file";
            return(-1);
        }

        //Push out value for this loop into the data vector. 
        t.push_back(std::stod(xVal));
        y.push_back(std::stod(yVal));
    }
    
    /*
    Check to make sure we have time series data. 
    */
   if (t.size() != y.size())
   {
       std::cout << "Error --- input data of different lengths";
       return(-1);
   }
   

    //Fit model
    Spline model(method, power, knots);
    model.fit(t,y);

    //Temp vectors to hold predictions. 
    std::vector<double> xTemp = linspace(t[0],t.back(), 0.1); //Model output.
    std::vector<double> modelTemp = {}; //Model output.

    //Predict values at each timepoint in our initial dataset. 
    for (uint64_t i = 0; i < xTemp.size(); i++)
    {
        modelTemp.push_back(model.predict(xTemp[i]));
    }
    
    /* 
    file output
    */
    std::ofstream output("output.csv");  //IO check.
    if (!output.is_open())
    {
        std::cout << "Error opening output file";
        return( -1);
    } 

    //Output CSV header names.
    output << "t" << ',';
    output << "y";
    output << '\n';

    //Populate .CSV.

    for (uint64_t i = 0; i < modelTemp.size(); i++)
    {
        output << xTemp[i] << ",";
        output << modelTemp[i] <<  '\n';
    }
    
    output.close();
    std::cout << "output.csv updated \n"; //

    PrintMat(model.Coe);

    std::vector<double> k2 = {0, 0, 0, 0, 2.8, 5.6, 8.4, 11.2, 14, 15, 15, 15, 15};

    std::cout << model.predict(8) << "\n";



    std::cout << "\n";
    std::cout << "done :D" << "\n";
    return 0;
}

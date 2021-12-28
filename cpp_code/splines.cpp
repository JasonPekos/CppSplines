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

    if(InputCheck(argc, argv) == 0){
        return(-1);
    }

    std::string method(argv[1]);

    if (InputCheck(argc, argv) == 2)
    {
        method = "PolynomialRegression"; //If we don't have any knot points, just use PNR as the default. 
    }
    std::string Arg2(argv[2]);
    std::string Arg3(argv[3]);

    uint64_t power = (uint64_t)stoi(Arg2);
    uint64_t knots = (uint64_t)stoi(Arg3);
    

    std::vector<double> t = {}; //Empty vectors to push back into
    std::vector<double> y = {};

    std::string xVal; //temp variables to store IOStream values.
    std::string yVal;

    /*
    Read in input.csv file
    */

    ifstream data("input.csv");
    if (!data.is_open()){
        cout << "Error opening file!";
        return(-1);
    }
    while (data.peek()!= EOF)
    {

        std::getline(data, xVal, ',');
        std::getline(data, yVal, '\n');


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


    std::vector<std::vector<double>> qq = model.Coe;


    /* 
    file output
    */
    std::vector<double> xTemp = t; //Model output.
    std::vector<double> modelTemp = {}; //Model output.

    for (uint64_t i = 0; i < t.size(); i++)
    {
        modelTemp.push_back(model.predict(t[i]));
    }
    


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


    for (uint64_t i = 0; i < modelTemp.size(); i++)
    {
        output << xTemp[i] << ",";
        output << modelTemp[i] <<  '\n';
    }
    
    output.close();
    std::cout << "output.csv updated \n"; //

    



    /*DO TESTING HERE*/

    Spline modelPowerBasis("PowerBasis", 3, 4);


    modelPowerBasis.fit(t,y);


    std::cout << "done :D" << "\n";
    return 0;
}

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

    //Bring over knots and power CL arguments.


    std::string PNR  = "PolynomialRegression";
    std::string BSP  = "BSpline";
    std::string PSP  = "PowerBasis";
    std::string SMS  = "Smooth";
    

    uint64_t knots = 0;
    uint64_t power = 3;
    double   lambda = 2;
    bool autoFlag = 0;
    //Add parameters from command line arguments. 
    if (method == SMS)
    {
        std::string Auto  = "auto";
        std::string Arg2(argv[2]);

        if (Arg2 == Auto)
        {
            autoFlag = 1;
        }
        else
        {
            lambda = stod(Arg2);
        }
    }
    
    if (method != SMS)
    {
        std::string Arg2(argv[2]);
        power = (uint64_t)stoi(Arg2);
    }

    if (method != PNR && method != SMS)
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
    //Skip header. 
    std::getline(data,xVal, ',');
    std::getline(data,yVal, '\n');
    if (xVal.find_first_not_of("-1234567890.") == std::string::npos)
    {
        std::cout << "warning! check to make sure .csv has column names. " << "\n";
    }
    if (yVal.find_first_not_of("-1234567890.") == std::string::npos)
    {
        std::cout << "warning! check to make sure .csv has column names. " << "\n";
    }
    
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
 
    //Temp vectors to hold predictions. 
    std::vector<double> xTemp = linspace(t[0],t.back(), 0.1); //Model output.
    std::vector<double> modelTemp = {}; //Model output.

    /* 
    FILE OUTPUT
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



    if(method != SMS)
    {
        Spline model(method, power, knots);

        model.fit(t,y);
        for (uint64_t i = 0; i < xTemp.size(); i++)
        {
            modelTemp.push_back(model.predict(xTemp[i]));
        }
    
    }
    else
    {
        GAM model(lambda);

        if (autoFlag == 1)
        {
            model.seekLambda(t,y);
            std::cout << "Recovered Lambda: "<< model.Lambda << "\n";
        }

        model.fit(t,y);
        for (uint64_t i = 0; i < xTemp.size(); i++)
        {
            modelTemp.push_back(model.predict(xTemp[i]));
        }
    }

    //Populate .CSV.
    for (uint64_t i = 0; i < modelTemp.size(); i++)
    {
        output << xTemp[i] << ",";
        output << modelTemp[i] <<  '\n';
    }
    
    output.close();
    std::cout << "output.csv updated \n"; //

    std::vector<double> kTemp = {4,8,12};

  

    std::vector<std::vector<double>> X = Design(t, 1);
    std::vector<std::vector<double>> XTX = MatMul(Transpose(X), X);
    std::vector<std::vector<double>> invXTX = Inverse(XTX);
    PrintMat(XTX);
    PrintMat(MatMul(invXTX, XTX));
    


    //PrintMat(model.Coe); Uncomment to return model coefficient values. For testing / debugging. 
    return 0;
}

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


    std::vector<double> told = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    std::vector<double> yold = {1,4,21,18,17,22,15,14,10,9,11,8,13,14,17};

    std::vector<double> t = {}; //Empty vectors to push back into
    std::vector<double> y = {}; 

    std::string line = "";
    /*
    Read in input.csv file
    */

    ifstream data("input.csv");
    if (!data.is_open()){
        cout << "Error opening file!";
        return(-1);
    }
    while (std::getline(data, line))
    {
        if (line.find_first_of("0123456789") == std::string::npos){
            continue;
        }

    }
    





    /*
    Check to make sure we have time series data. 

    */
    


    
    





    // std::string method(argv[1]);

    // Spline model("PolynomialRegression",3, 0);

    // model.fit(t,y);
    // double a = model.predict(12);

    // std::cout << "predicted value:" << a << "\n";

    // std::vector<std::vector<double>> qq = model.Coe;

    // PrintMat(qq);



    /* 
    file output
    */
    std::vector<double> xTemp = {1,4,21,18,17,22,15,14,10,9,11,8,13,14,17}; //Model output.
    std::vector<double> modelTemp = {1,4,21,18,17,22,15,14,10,9,11,8,13,14,17}; //Model output.


    std::ofstream output("output.csv");  //IO check.
    if (!output.is_open())
    {
        std::cout << "Error opening output file";
        return( -1);
    } 

    //Output CSV header names.
    output << "t" << ',';
    output << "y" << ',';
    output << '\n';


    for (uint64_t i = 0; i < modelTemp.size(); i++)
    {
        output << xTemp[i] << ",";
        output << modelTemp[i] << "," <<  '\n';
    }
    
    output.close();
    std::cout << "output.csv updated \n"; //

    std::cout << "done :D" << "\n";

    return 0;
}

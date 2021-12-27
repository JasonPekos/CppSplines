#include "model.hpp" //Leave out until LA is working
#include "linAlgHelperFunctions.hpp" // when I move everything over.
#include "iocheck.hpp" //iochecking
#include <vector>
#include <string>
#include <iostream>
#include <cmath>



int main(int argc, char const *argv[])
{

    if(InputCheck(argc, argv) == 0){
        return(-1);
    }


    
    

    std::vector<double> t = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    std::vector<double> y = {1,4,21,18,17,22,15,14,10,9,11,8,13,14,17};



    std::string method(argv[1]);

    Spline model("PolynomialRegression",3, 0);

    model.fit(t,y);
    double a = model.predict(12);

    std::cout << "predicted value:" << a << "\n";

    std::vector<std::vector<double>> qq = model.Coe;

    PrintMat(qq);


    std::cout << "done :D" << "\n";

    return 0;
}

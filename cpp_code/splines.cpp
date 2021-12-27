//#include "model.hpp" Leave out until LA is working
#include "linAlgHelperFunctions.hpp" // when I move everything over.
#include <vector>
#include <string>
#include <iostream>
#include <cmath>



int main(int argc, char const *argv[])
{
    if (argc < 2)
    {
        std::cout << "add method" << "\n";
        exit(-1);
    }
    


    std::string method(argv[1]);
    std::string PNR("PolynomialRegression");

    std::vector<double> t = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    std::vector<double> y = {1,4,21,18,17,22,15,14,10,9,11,8,13,14,17};



    if (method == PNR)
    {
        std::cout << method << "\n";
        std::vector<std::vector<double>> mat = Design(t, 4);

        std::vector<std::vector<double>>  matT = Transpose(mat);

        std::vector<std::vector<double>> Soln = MatVecMul(matT, y);
        //PrintMat(GaussElimination(matT,y1));

        PrintMat(Soln);

        //TriLCheck(matT);
        std::cout << "Design Matrix:" << "\n";
        //PrintMat(mat);

        std::cout << "Q'Q:" << "\n";
        //PrintMat(MatMul(matT,mat));

        std::vector<std::vector<double>> Z = SolveSystem(MatMul(matT,mat), Soln);

        PrintMat(Z);

        PrintMat(MatMul(matT, mat));

        // coe = (X' * X ) \ (X'y)
         //Need: matrix transpose, matrix multiplication, solve linear system 
    }
    


    std::vector<std::vector<double>> MatrixA(20, std::vector<double>(20,0));
    std::vector<std::vector<double>> MatrixB(20, std::vector<double>(20,0));

    /* data IO */

    //if iocheck() == 1, cont, otherwise exit


    /* actual model usage*/

    //Instantiate 
    // Spline sModel(basis, power, knots)

    
    //specific_spline.fit(data)

    //specific_spline.predict(x)

    std::cout << "done :D" << "\n";

    return 0;
}

#include<iostream>   // Error messages. 
#include<string>     //String processing.
#pragma once



int InputCheck(int argc, const char** argv){

    /**
     * @brief This function checks all the inputs except the validity of the methods used.
     * 
     * @param argc number of arguments, standard from C++ main() call when using input at command line.
     * 
     * @param argv arguments character pointer, standard from C++ main() call when using input at command lin
     * 
     */
    
    if (argc < 3) // Checks to make sure we have at least three arguments.
    {
        std::cout << "There is an issue here with the number of parameters submitted — ";
        std::cout << "the minimum number of parameters for any procedure is 2 \n";
        std::cout << "\n Parameter One: Method (e.g PolynomialRegression) \n";
        std::cout << " Parameter Two: degree of spline basis (e.g. 3)\n";
        std::cout << " Parameter Three (for BSpline and PowerBasis Spline): knots (e.g. 2) \n";
        return(0);
    }

    //Read in method string because we need to check against method to see what we need to look at next.
    std::string method(argv[1]);
    std::string PNR  = "PolynomialRegression";
    std::string BSP  = "BSpline";
    std::string PSP  = "PowerBasis";


    /*
    Check if input parameters are numbers where they should be 
    (don't need to check the converse because making into string and checking against valid inputs handles that)
    */
    std::string Arg2(argv[2]);

    if (method != PNR)
    {
        if (argc < 4) // Checks to make sure we have at least three arguments.
        {
            std::cout << "There is an issue here with the number of parameters submitted — ";
            std::cout << "the minimum number of parameters for this procedure is 3 \n";
            std::cout << "\n Parameter One: Method (e.g BSpline) \n";
            std::cout << " Parameter Two: degree of spline basis (e.g. 3)\n";
            std::cout << " Parameter Three (for BSpline and PowerBasis Spline): knots (e.g. 2) \n";
            return(0);
        }
        std::string Arg3(argv[3]);
        if (Arg3.find_first_not_of("1234567890") != std::string::npos)
        {
            std::cout << "Error: Please submit a positive integer as argument three.\n";
            return(0);
        }
    }
    
    if (Arg2.find_first_not_of("1234567890") != std::string::npos) //Check to make sure arg two (power) is a positive integer. 
    {
        std::cout << "Error: Please submit a positive integer as argument two.\n";
        return(0);
    }

    
    //Check  arguments for possible errors and warnings.
    if (method != PNR)
    {
        if (method != PSP)
        {
            if (method != BSP)
            {
                std::cout << "Please submit a valid method (e.g. PolynomialRegression)" << "\n";
                exit(-1);
            } 
        }
    }
    else if (atof(argv[2]) > 100) //Splines of degree this large are very dumb. 
    {
        std::cout << "You have submitted an unreasonably large power value. ";
        std::cout << "The code will run, and this is mathematically well defined but please ";
        std::cout << "consider trying a smaller value --- splines of degree three are standard! \n";
        //return(0);
    }
    else if (method != PNR) //Check number of knots submitted. 
    {
        if (atof(argv[3]) == 0)
        {
            std::cout << "Zero knots submitted.\n";
            std::cout << "Changing to polynomial regression. \n";
            return(2);
        }
    }
    return(1);
}
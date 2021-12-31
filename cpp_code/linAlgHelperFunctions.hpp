#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "mathHelperFunctions.hpp" // when I move everything over.

#pragma once

/*
This file contains linear algebra helper functions to help with various computations:

1. PrintMat: prints a given matrix (vector of vectors) to the terminal.
2. Design: Creates a Vandermone matrix 
3. Transpose: Return matrix transpose
4. MatMul: Matrix multiplication
5. MatVecMul: Matrix-Vector multiplication + standardizing type.
6. TriLCheck: Check if a matrix is lower triangular. 
7. SolveSystem: Solve a linear system using reduction to Row Echelon form and then back substitution.
*/


void PrintMat(std::vector<std::vector<double>> A){
    /**
     * @brief This function prints a matrix to the terminal for troubleshooting.
     * 
     * @param A Matrix to be printed. 
     * 
     */

    //This is here to break up different print calls.
    std::cout << "requested matrix:" << "\n";  

    //For each row, print every column, and then advance. 
    for (uint64_t i = 0; i < A.size(); i++) 
    {
        for (uint64_t j = 0; j < A[1].size(); j++)
        {
            std::cout << A[i][j] << " "; 
        }   

    //New line for formatting.
    std::cout << "\n";  
    }
}

std::vector<double> StripDuplicates(std::vector<double> a){
    /**
     * @brief Strip duplicates that are next to eachother from a vector.
     * 
     * @param a Input vector
     * 
     * @return Input vector without consecutive duplicates. 
     * 
     */
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    std::vector<double> out = a;
    return(out);
}

std::vector<std::vector<double>> Design(std::vector<double> t, uint64_t power){
    /**
     * @brief Creates a Vandermonde matrix:
     * 
     * https://en.wikipedia.org/wiki/Vandermonde_matrix
     * 
     * For fitting polynomial regression models. 
     * 
     * @param t Time vector from data of form (t,y)
     * @param power Largest element has degree power - 1.
     * 
     * @return length(t) by j Vandermonde corresponding to input data.  
     */

    power = power + 1;
    //Set up matrix with correct dimensions:
    std::vector<std::vector<double>> mat(t.size(), std::vector<double>(power)); 

    //Populate matrix with appropriate elements: mat_i,j = t_i^{j-1}
    for (uint64_t i = 0; i < t.size(); i++)
    {
        for (uint64_t j = 0; j < power; j++)
        {
            mat[i][j] = pow(t[i],j);
        }
        
    }
    return(mat);
}

std::vector<std::vector<double>> DesignPowerBasis(std::vector<double> t, uint64_t power, std::vector<double> knots){
    /**
     * @brief Creates a design matrix for fitting a power basis spline of the nth order:
     * 
     * 
     * For fitting polynomial regression models. 
     * 
     * @param t Time vector from data of form (t,y)
     * @param power Largest element has degree power - 1.
     * @param knots Number of knots
     * 
     * @return length(t) by j + k corresponding to input data.  
     */

    power = power + 1;

    //Set up matrix with correct dimensions:
    std::vector<std::vector<double>> mat(t.size(), std::vector<double>(power + knots.size())); 

    //Populate matrix with appropriate elements: mat_i,j = t_i^{j-1}
    for (uint64_t i = 0; i < t.size(); i++)
    {
        for (uint64_t j = 0; j < power; j++)
        {
            mat[i][j] = pow(t[i],j);
        }
        for (uint64_t j = 0; j < knots.size(); j++) //
        {
            mat[i][power + j] = pow(pm(t[i] - knots[j]),(power - 1));
        }
    }
    return(mat);
}

std::vector<std::vector<double>> DesignBSplineBasis(std::vector<double> t, uint64_t power, std::vector<double> knots){
    /**
     * @brief Creates a design matrix for fitting a B-spline of the nth order:
     * 
     * 
     * For fitting regression spline models. 
     * 
     * @param t Time vector from data of form (t,y)
     * @param power Largest element has degree power - 1.
     * @param knots Vector of knot positions. 
     * @param InteriorKnots Number of interior knots
     * 
     * @return length(t) by j + k corresponding to input data.  
     */

    //power = power + 1;
    std::vector<double> RealKnots = StripDuplicates(knots);
    uint64_t k = RealKnots.size() + power + 1 - 2;

    //Set up matrix with correct dimensions:
    std::vector<std::vector<double>> mat(t.size(), std::vector<double>(k, 0));

    
    //Populate matrix with appropriate elements: 
    for (uint64_t i = 0; i < t.size(); i++)
    {
        for (uint64_t j = 1; j <= k; j++)
        {
            mat[i][j - 1] = CoxDeBoor(t[i], j , knots, power);   
        }
    }
    mat[t.size() - 1][k-1] = 1; //Assign at right boundary to properly fit boundary basis. 
    return(mat);
}


std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> A){

    /**
     * @brief Transpose matrix
     * 
     * @param A Matrix to be transposed
     * 
     * @return Tranposed version of A
     */


    std::vector<std::vector<double>> B(A[1].size(), std::vector<double>(A.size())); //Create new matrix with identical dimensions to A.

    //Reassign values. 
    for (uint64_t i = 0; i < A.size(); i++)
    {
        for (uint64_t j = 0; j < A[1].size(); j++)
        {
            B[j][i] = A[i][j];
        }
    }
    //Return _new_ matrix. 
    return(B);
}

std::vector<std::vector<double>> MatMul(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){
    /**
     * @brief Multiply two matrices together and return a new matrix. 
     * 
     * @param A First matrix
     * @param B Second matrix
     * 
     * @return C = A*B.
     */

    //Create new matrix of size A_i, B_j, populate with zeroes. 
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(B[1].size(), 0));

    //Loop over all A rows and B columns and add the dot product to the output matrix.
    for (uint64_t i = 0; i < A.size(); i++)
    {
        for (uint64_t j = 0; j < B[1].size(); j++)
        {
            double sum = 0;
            for (uint64_t k = 0; k < A[1].size(); k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;    
        } 
    }
    //Return output as new matrix. 
    return(C);
}

std::vector<std::vector<double>> MatVecMul(std::vector<std::vector<double>> A, std::vector<double> y){

    /**
     * @brief Multiply a matrix by a vector for the right.  
     * 
     * @param A Input matrix in Ay = x
     * @param y Input vector in Ay = x
     * 
     * @return x as a _matrix_, not a vector, to keep classes consistent in later use. E.g. [y] instead of y. 
     */

    //Create y, populate with zeroes. 
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(1, 0));

    for (uint64_t i = 0; i < A.size(); i++)
    {

        double sum = 0;
        for (uint64_t k = 0; k < y.size(); k++)
        {
            sum += A[i][k] * y[k]; //Push vector through matrix and add result to output. 
        }
        C[i][0] = sum;    
    }
    //Return C = Ay. 
    return(C);
}

bool TriLCheck(std::vector<std::vector<double>> A){
    /**
     * @brief Check if a matrix is upper triangular. E.g. all element below the primary diagonal are zero. 
     * 
     * @param A Matrix in question
     * 
     * @return Bool: one if matrix is Upper Triangular, zero else.
     */

    double sum = 0;
    for (uint64_t i = 1; i < A.size(); i++) //start at one because we don't care about first row. 
    {
        for (uint64_t j = 0; j < i; j++) // Sum absolute value of all elements to see if any are greater than zero. 
        {
            sum += abs(A[i][j]);
        }   
    }
    if (sum == 0 ) //Matrix is Lower Triangular 
    {
        return(1);
    }
    else //There is a nonzero element below the main diagonal.
    {
        return(0);
    }
}

std::vector<std::vector<double>> SolveSystem(std::vector<std::vector<double>> B, std::vector<std::vector<double>> y){

    /**
     * @brief Solve a linear system: find 'a' such that Ba = y. 
     * 
     * @param B Corresponding matrix
     * @param y Given data, e.g. y in (t,y) time series. 
     * 
     * @return 1d Matrix of coefficients that solve the linear system; 'a' in Ba = y. 
     */

    //Create new vectors of the size that we want by copying the input values. 
    std::vector<std::vector<double>> V = y; //For use augmenting matrix. 
    std::vector<std::vector<double>> Soln = y; //For use in Back Sub
    std::vector<std::vector<double>> A = B; //For use in both algorithms. 

    ///Gaussian Elimination

    //Check if we have a unique solution --- if not, a transpose is missing somewhere, and we need to throw an error to avoid segfaulting. 
    if (A.size() < A[1].size())
    {
        std::cout << "Overdetermined System of Linear Equations --- Error.";
        exit(-1);
    }

    //Define values for the While loop
    uint64_t maxiters = 1000000; //Really big number to avoid a totally fatal infinite loop. 
    uint64_t count    = 0; //Count for while loop.

    //Augment Matrix.
    for (uint64_t i = 0; i < A.size(); i++)
    {
        A[i].push_back(V[i][0]);
    }

    //Start actual Gaussian elimination process. 
    while (count < maxiters)
    {
        for (uint64_t index = 0; index < A.size(); index++)
        {
            //Loop over all the indices from zero -> number of rows (not columns; matrix is augmented and no longer square)
            for (uint64_t rowIndex = index + 1; rowIndex < A.size(); rowIndex++){
                if (A[rowIndex][index] != 0) //If already zero we _must_ skip to avoid division by zero.
                {
                    double currentVal = A[rowIndex][index]; //Values for computation after quick division by zero check. 
                    double rowOneVal = A[index][index]; 

                    if (rowOneVal != 0)
                    {
                        double correctionFactor = -(currentVal / rowOneVal);
                        //Add a multiple of correctionfactor*rowoneindex + rowindex DO THIS
                        for (uint64_t iter = 0; iter < A[1].size(); iter++)
                        {
                            A[rowIndex][iter] = A[rowIndex][iter] + correctionFactor*A[index][iter]; //Add row*multiple to row to reduce to zero. 
                        }
                    }
                }
            }
        }
        if (TriLCheck(A) == 1) //Check if our matrix is lower triangular -> simplified. If so, break. 
        {
            count = maxiters;
        }
        count = count + 1;        
    }


    
    //Back Substitution.

    for (int64_t i =  (int64_t)A.size() - 1; i >= 0; --i) //Starting at last element, decrementing. In back sub we roll up from the last (simplist) row.
    {
        uint64_t m_size = A[1].size(); //Define outside of indexing into an array for readability.

        double current = A[(uint64_t)i][m_size - 1]; //Think about last row of a lower triangular matrix; we want to start with [0,0,0 .., 0, x | y], so grab len - 1th element. 

        double cumsum = 0;
        for (uint64_t j = (uint64_t)i+1; j < A.size(); j++)
        {
            cumsum = cumsum + A[(uint64_t)i][j] *  Soln[(uint64_t)j][0]; //See e.g. https://algowiki-project.org/en/Backward_substitution 1.2
        }

        double answer =  ((current - cumsum) / A[(uint64_t)i][(uint64_t)i]);
        Soln[(uint64_t)i][0] = answer;
    }
    return(Soln); //Return solution vector.
}

bool NotNAN(double var){
    /**
     * @brief Checks if a value is NAN
     * 
     * @param var Variable we're checking.
     * 
     */
    if (var != var)
    {
        return(0);
    }
    else
    {
        return(1);
    }
}

bool MatNoNAN(std::vector<std::vector<double>> mat){
    /**
     * @brief Returns One if Vector contains a NaN.
     * 
     * @param mat Matrix we care about. 
     * 
     */
    for (uint64_t i = 0; i < mat.size(); i++)
    {
        for (uint64_t j = 0; j < mat[0].size(); j++)
        {
            if (!NotNAN(mat[i][j]))
            {
                return(0);
            }
        }
    }
    return(1);
}

std::vector<std::vector<double>> Eye(uint64_t size1, uint64_t size2){
    /**
     * @brief Returns an identity matrix of size 'size1, size2'. Useful for constructing penalization matrices. 
     * 
     * @param size Function returns a [size1, size2] identity matrix.
     * 
     * @return size1 by size2 identity matrix.
     * 
     */

    //Create output, populate with zeros.
    std::vector<std::vector<double>> out(size1, std::vector<double>(size2, 0));

    //Add ones along diagonal.
    for (uint64_t i = 0; i < out.size(); i++)
    {
        for (uint64_t j = 0; j < out[0].size(); j++)
        {
            if (i == j )
            {
                out[i][j] = 1;
            }
        }       
    }
    return(out);
}


std::vector<std::vector<double>> AddWigglyPenalty(double lambda, std::vector<std::vector<double>> B){
    /**
     * @brief This function adds a wigglyness penalty to a smoothing spline. The derivation
     * for this penalty matrix is given in Wood — Generalized Additive Models (2006).
     * 
     * It provides a slightly crude approximation to the integral of the second derivative.
     * 
     * Issues can arise when data is highly irregularly spaced out, though this was already an issue with smooths,
     * given automatic knot location.  
     * 
     * @param lambda This is the wigglyness parameter for the model, determining how much we let the penalization matrix flatten
     * our fit. As lambda -> infinity the model becomes linear; as lambda -> 0, the model becomes a perfect cubic spline interpolation. 
     * 
     * @param B This is the matrix --- usually XTX in a regression problem --- that we want to penalize. 
     * 
     * @return A matrix with dimensions equivalent to B, with the penalty added for a generalized Ridge Regression problem of type:
     * BetaHat = (XTX + penalty) \ XT y.  
     * 
     */
    
    std::vector<std::vector<double>> out = B;

    //Initialize a vector of zeroes at the correct size (dimensions given in Wood 2006).
    std::vector<std::vector<double>> D = Eye(B.size() - 2, B.size());

    for (uint64_t i = 0; i < D.size(); i++)
    {
        for (uint64_t j = 0; j < D[0].size(); j++)
        {
            if (i == j)
            {
                //Add one along the diagonal.
                D[i][j] = 1;
            }
            if (i + 1 == j)
            {
                //Add -2 on the first off diagonal.
                D[i][j] = -2;
            }
            if (i + 2 == j)
            {
                //Add 1 on the second off diagonal.
                D[i][j] = 1;
            }
        }
    }

    //Matrix multiplication to recover S = DTD.
    D = MatMul(Transpose(D), D);

    //Add penalty to input matrix. 
    for (uint64_t i = 0; i < out.size(); i++)
    {
        for (uint64_t j = 0; j < out.size(); j++)
        {
            out[i][j] += pow(lambda,2) * D[i][j];
        }
    }

    //Returns (Input + lambda^2 DTD).
    return(out);
}

bool IsEye(std::vector<std::vector<double>> A, uint64_t l1, uint64_t l2){
    /**
     * @brief Checks if a matrix is an identity matrix. l1, l2 used to check section of an augmented matrix. 
     * 
     * Normal use: IsEye(A, A.size(), A[0].size())
     * 
     * @param A Matrix to be checked.
     * @param l1 Number of rows.
     * @param l2 Number of columns. 
     * 
     * @return True if matrix is an identity
     * 
     */
    for (uint64_t i = 0; i < l1; i++)
    {
        for (uint64_t j = 0; j  < l2; j ++)
        {
            if (i != j)
            {
                if (A[i][j] != 0)
                {
                    return(0);
                }
            }
            if (i == j)
            {
                if (A[i][j] != 1)
                {
                    return(0);
                }     
            }
        }
    }
    return(1);
}


std::vector<std::vector<double>> Inverse(std::vector<std::vector<double>> A){
    /**
     * @brief Calculates the matrix inverse. 
     * 
     */

    if (A.size() != A[0].size())
    {
        std::cout << "error, matrix is not square!" << "\n";
    }
    

    std::vector<std::vector<double>> IDtoAugment = Eye(A.size(), A[0].size());
    std::vector<std::vector<double>> Inv(A.size(), std::vector<double>((A[0].size() + IDtoAugment[0].size()), 0));
    std::vector<std::vector<double>> C = A;


    //Augment Matrix.
    for (uint64_t i = 0; i < A.size(); i++)
    {
        for (uint64_t j = 0; j < A[0].size(); j++)
        {
            Inv[i][j] = A[i][j];
        }
        for (uint64_t k = A[0].size(); k < (A[0].size() + IDtoAugment[0].size()); k++)
        {
            Inv[i][k] = IDtoAugment[i][k - A[0].size()];
        }
    }
    //Augmentation Complete.

     //Define values for the While loop
    uint64_t maxiters = 1000000; //Really big number to avoid a totally fatal infinite loop. 
    uint64_t count    = 0; //Count for while loop.


    //Start actual Gaussian elimination process. 
    while (count < maxiters)
    {
        for (uint64_t index = 0; index < Inv[0].size(); index++)
        {
            //Loop over all the indices from zero -> number of rows (not columns; matrix is augmented and no longer square)
            for (uint64_t rowIndex = 0; rowIndex < A.size(); rowIndex++){

                //Inv[rowIndex][index] = Inv[rowIndex][index] / Inv[rowIndex][rowIndex];

                if (Inv[rowIndex][index] != 0) //If already zero we _must_ skip to avoid division by zero.
                {
                    double currentVal = Inv[rowIndex][index]; //Values for computation after quick division by zero check. 
                    double rowOneVal = Inv[index][index]; 

                    if (rowOneVal != 0)
                    {
                        if (rowIndex != index)
                        {
                            double correctionFactor = -(currentVal / rowOneVal);
                            //Add a multiple of correctionfactor*rowoneindex + rowindex DO THIS
                            for (uint64_t iter = 0; iter < Inv[0].size(); iter++)
                            {
                                PrintMat(Inv);
                                Inv[rowIndex][iter] = Inv[rowIndex][iter] + correctionFactor*Inv[index][iter]; //Add row*multiple to row to reduce to zero. 
                                
                            }
                        }
                    }
                }
            }
        }

        //Make principle diagonal into one. 
        for (uint64_t i = 0; i < Inv.size(); i++)
        {
            for (uint64_t j = 0; j < Inv[0].size(); j++)
            {
                
            }   
        }
        
  
        //Check if our matrix is Identity. If so, break. 
        if (IsEye(Inv, A.size(), A[0].size()) == 1) 
        {
            count = maxiters;
        }
        count = count + 1;        
    }

    //Check if LHS is identity:

    std::vector<std::vector<double>> Out = A;

    for (uint64_t i = A.size(); i < Inv.size(); i++)
    {
        for (uint64_t j = A[0].size(); j < Inv[0].size(); j++)
        {
            Out[i][j] = Inv[i][j];
        } 
    }
    
    return(Out);
}

double Rank(std::vector<std::vector<double>> A){
    /**
     * @brief Returns the rank for some matrix A.
     * 
     */

    double sum = 0;
    for (uint64_t i = 0; i < A.size(); i++)
    {
        for (uint64_t j = 0; j < A[0].size(); j++)
        {
            if (i == j)
            {
                 sum += A[i][j];
            }
        }
    }
    return(sum);
}
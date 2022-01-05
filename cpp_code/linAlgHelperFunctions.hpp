/**
 * @file linAlgHelperFunctions.hpp
 * @author Jason Pekos
 * @brief This file contains linear algebra helper functions to help with various computations:
 * 
 * 1.  PrintMat: prints a given matrix (vector of vectors) to the terminal.
 * 2.  Design: Creates a Vandermone matrix 
 * 3.  Transpose: Return matrix transpose
 * 4.  MatMul: Matrix multiplication
 * 5.  MatVecMul: Matrix-Vector multiplication + standardizing type.
 * 6.  TriUCheck: Check if a matrix is upper triangular. 
 * 7.  SolveSystem: Solve a linear system using reduction to Row Echelon form and then back substitution.
 * 8.  StripDuplicates: Strip duplicates that are next to eachother from a vector.
 * 9.  DesignPowerBasis: Vandermonde style design matrix for a power basis problem.
 * 10. DesignBSplineBasis: Vandermonde style design matrix for a BSpline problem.
 * 11. NotNAN / MatNotNAN: NAN checking for use in regression splines.
 * 12. AddWigglyPenalty: Add the penalty --- lambda^2 \int f''(x) dx --- to a regression problem.
 * 14. Eye: Create identity matrix.
 * 13. IsEye: Check if a matrix is the identity matrix.
 * 14. Inverse: Invert a matrix.
 * 15. Trace: Find the trace of a matrix.
 * 16. TraceHatMatrix: Find the trace of a hat matrix in the context of a GAM problem.
 * 
 * @version 0.1
 * @date 2021-12-31
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <vector>                   // Everything here uses the vector class. 
#include <string>                   // Strings for holding methods etc.
#include <iostream>                 // Debugging, printing to terminal, check-ins for long algorithms. 
#include <cmath>                    // Pow() etc.
#include <algorithm>                // Iterators for vector sort. 
#include "mathHelperFunctions.hpp"  // Cox-De-Boor construction, pm for power basis.

#pragma once


/**
 * @brief This function prints a matrix to the terminal for troubleshooting.
 * 
 * @param A Matrix to be printed. 
 */
void PrintMat(std::vector<std::vector<double>> A)
{
    // This is here to break up different print calls.
    std::cout << "requested matrix:"
              << "\n";

    // For each row, print every column, and then advance.
    for (uint64_t i = 0; i < A.size(); i++)
    {
        for (uint64_t j = 0; j < A[1].size(); j++)
        {
            std::cout << A[i][j] << " ";
        }

        // New line for formatting.
        std::cout << "\n";
    }
}


/**
 * @brief Strip duplicates that are next to eachother from a vector.
 * 
 * @param a Input vector
 *
 * @return Input vector without consecutive duplicates. 
 */
std::vector<double> StripDuplicates(std::vector<double> a)
{
    // Call sort with iterators. we are using this for knot vectors, so they ought to be sorted anyways. 
    std::sort(a.begin(), a.end());

    // Strip duplicates. 
    a.erase(std::unique(a.begin(), a.end()), a.end());
    std::vector<double> out = a;
    return (out);
}


/**
 * @brief Creates a Vandermonde matrix:
 * 
 * https:// en.wikipedia.org/wiki/Vandermonde_matrix
 * 
 * For fitting polynomial regression models. 
 * 
 * @param t Time vector from data of form (t,y)
 * @param power Largest element has degree power - 1.
 * 
 * @return length(t) by j Vandermonde corresponding to input data.  
 */
std::vector<std::vector<double>> Design(std::vector<double> t, uint64_t power)
{
    power = power + 1;
    // Set up matrix with correct dimensions:
    std::vector<std::vector<double>> mat(t.size(), std::vector<double>(power));

    // Populate matrix with appropriate elements: mat_i,j = t_i^{j-1}
    for (uint64_t i = 0; i < t.size(); i++)
    {
        for (uint64_t j = 0; j < power; j++)
        {
            mat[i][j] = pow(t[i], (double)j);
        }
    }
    return (mat);
}


/**
 * @brief Creates a design matrix for fitting a power basis spline of the nth order.
 * Used primarily for fitting power basis regression models. 
 * 
 * @param t Time vector from data of form (t,y)
 * @param power Largest element has degree power - 1.
 * @param knots Number of knots
 * 
 * @return length(t) by j + k corresponding to input data.  
 */
std::vector<std::vector<double>> DesignPowerBasis(std::vector<double> t, uint64_t power, std::vector<double> knots)
{
    // Power n := order n+1 spline.
    power = power + 1;

    // Set up matrix with correct dimensions:
    std::vector<std::vector<double>> mat(t.size(), std::vector<double>(power + knots.size()));

    // Populate matrix with appropriate elements: mat_i,j = t_i^{j-1}
    for (uint64_t i = 0; i < t.size(); i++)
    {
        for (uint64_t j = 0; j < power; j++)
        {
            mat[i][j] = pow(t[i], (double)j);
        }
        for (uint64_t j = 0; j < knots.size(); j++) // 
        {
            mat[i][power + j] = pow(pm(t[i] - knots[j]), ((double)power - 1));
        }
    }
    return (mat);
}


/**
 * @brief Creates a design matrix for fitting a B-spline of the nth order:
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
std::vector<std::vector<double>> DesignBSplineBasis(std::vector<double> t, uint64_t power, std::vector<double> knots)
{
    /*
    This matrix is fed knots beyond the edge of the domain, due to the convolution based method for B-Spline construction
    requiring padded knots. Here we strip those away before creating the design matrix.
    */
    std::vector<double> RealKnots = StripDuplicates(knots);
    uint64_t k = RealKnots.size() + power + 1 - 2;

    // Set up matrix with correct dimensions:
    std::vector<std::vector<double>> mat(t.size(), std::vector<double>(k, 0));

    // Populate matrix with appropriate elements:
    for (uint64_t i = 0; i < t.size(); i++)
    {
        for (uint64_t j = 1; j <= k; j++)
        {
            mat[i][j - 1] = CoxDeBoor(t[i], j, knots, power);
        }
    }
    // Assign at right boundary to properly fit boundary basis.
    mat[t.size() - 1][k - 1] = 1; 
    return (mat);
}


/**
 * @brief Transpose matrix
 * 
 * @param A Matrix to be transposed
 * 
 * @return Tranposed version of A
 */
std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> A)
{
    std::vector<std::vector<double>> B(A[1].size(), std::vector<double>(A.size())); // Create new matrix with identical dimensions to A.

    // Reassign values.
    for (uint64_t i = 0; i < A.size(); i++)
    {
        for (uint64_t j = 0; j < A[1].size(); j++)
        {
            B[j][i] = A[i][j];
        }
    }
    // Return _new_ matrix.
    return (B);
}


/**
 * @brief Multiply two matrices together and return a new matrix. 
 * 
 * @param A First matrix
 * @param B Second matrix
 * 
 * @return C = A*B.
 */
std::vector<std::vector<double>> MatMul(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
    // Create new matrix of size A_i, B_j, populate with zeroes.
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(B[1].size(), 0));

    // Loop over all A rows and B columns and add the sum to the output matrix.
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
    // Return output as new matrix.
    return (C);
}


/**
 * @brief Multiply a matrix by a vector to the right.  
 * 
 * @param A Input matrix in Ay = x
 * @param y Input vector in Ay = x
 * 
 * @return x as a _matrix_ (vector of vectors, of size (n,0)), not a vector, to keep types consistent in later use.
 */
std::vector<std::vector<double>> MatVecMul(std::vector<std::vector<double>> A, std::vector<double> y)
{
    // Create output, populate with zeroes.
    std::vector<std::vector<double>> C(A.size(), std::vector<double>(1, 0));

    for (uint64_t i = 0; i < A.size(); i++)
    {

        // Push vector through matrix and add result to output.
        double sum = 0;
        for (uint64_t k = 0; k < y.size(); k++)
        {
            sum += A[i][k] * y[k]; 
        }
        // Assign  sum at each index to the output vector.
        C[i][0] = sum;
    }
    // Return C = Ay.
    return (C);
}


/**
 * @brief Check if a matrix is upper triangular. E.g. all element below the primary diagonal are zero. 
 * 
 * @param A Matrix in question
 * 
 * @return Bool: one if matrix is Upper Triangular, zero else.
 */
bool TriUCheck(std::vector<std::vector<double>> A)
{
    double sum = 0;

    // Start at one because we don't care about first row.
    for (uint64_t i = 1; i < A.size(); i++) 
    {
        // Sum absolute value of all elements to see if any are greater than zero.
        for (uint64_t j = 0; j < i; j++) 
        {
            double val = A[i][j];
            sum += pow(val,2);
        }
    }
    /*
    Use this tolerance here because we need to use an inverse later to calculate the Hat matrix, which is ill conditioned. 
    Ideally, I would have implemented some decomposition to avoid this. 

    Issues around numerical algebra broke the solver if we checked for exact equivalency, so I needed to check within a tolerance. 
    */
    if (sum < 1.0e-200) 
    {
        return (1);
    }
    else // There is a nonzero element below the main diagonal, return FALSE.
    {
        return (0);
    }
}



/**
 * @brief Solve a linear system: find 'a' such that Ba = y. 
 * 
 * Uses Gaussian Elimination with Back Substitution. 
 * 
 * @param B Corresponding matrix
 * @param y Given data, e.g. y in (t,y) time series. 
 * 
 * @return 1d Matrix of coefficients that solve the linear system; 'a' in Ba = y. 
 */
std::vector<std::vector<double>> SolveSystem(std::vector<std::vector<double>> B, std::vector<std::vector<double>> y)
{
    // Create new vectors of the size that we want by copying the input values.
    std::vector<std::vector<double>> V = y;    // For use augmenting matrix.
    std::vector<std::vector<double>> Soln = y; // For use in Back Sub
    std::vector<std::vector<double>> A = B;    // For use in both algorithms.

    // Gaussian Elimination

    // Check if we have a unique solution --- if not, a transpose is missing somewhere, and we need to throw an error to avoid segfaulting.
    if (A.size() < A[1].size())
    {
        std::cout << "Overdetermined System of Linear Equations --- Error.";
        exit(-1);
    }

    // Define values for the While loop
    uint64_t maxiters = 1000000; // Really big number to avoid a totally fatal infinite loop.
    uint64_t count = 0;          // Count for use in while loop.

    // Augment Matrix with response vector.
    for (uint64_t i = 0; i < A.size(); i++)
    {
        A[i].push_back(V[i][0]);
    }

    // Start actual Gaussian elimination process.
    while (count < maxiters)
    {
        for (uint64_t index = 0; index < A.size(); index++)
        {
            // Loop over all the indices from zero -> number of rows (not columns; matrix is augmented and no longer square)
            for (uint64_t rowIndex = index + 1; rowIndex < A.size(); rowIndex++)
            {
                // If already zero we skip to avoid changing by nothing.
                if (A[rowIndex][index] != 0) 
                {
                    if (A[index][index] != 0)
                    {
                        double correctionFactor = -(A[rowIndex][index] / A[index][index]);
                        // Add a multiple of correctionfactor*rowoneindex + rowindex
                        for (uint64_t iter = 0; iter < A[1].size(); iter++)
                        {
                            // Add row*multiple to row to reduce to zero.
                            A[rowIndex][iter] = A[rowIndex][iter] + correctionFactor * A[index][iter]; 
                        }
                    }
                }
            }
        }
        // Check if our matrix is lower triangular -> simplified. If so, break.
        if (TriUCheck(A) == 1) 
        {
            count = maxiters;
        }
        count = count + 1;
    }

    /*
    Back Substitution.
    */

   // Starting at last element, decrementing. In back sub we roll up from the last (simplest) row.
    for (int64_t i = (int64_t)A.size() - 1; i >= 0; --i) 
    {
        // Define outside of indexing into an array for readability.
        uint64_t m_size = A[1].size(); 

        // Think about last row of a lower triangular matrix; we want to start with [0,0,0 .., 0, x | y], so grab len - 1th element.
        double current = A[(uint64_t)i][m_size - 1]; 

        double cumsum = 0;
        for (uint64_t j = (uint64_t)i + 1; j < A.size(); j++)
        {
            cumsum = cumsum + A[(uint64_t)i][j] * Soln[(uint64_t)j][0]; // See e.g. https:// algowiki-project.org/en/Backward_substitution 1.2
        }

        // Push answer at each step into solution vector
        double answer = ((current - cumsum) / A[(uint64_t)i][(uint64_t)i]);
        Soln[(uint64_t)i][0] = answer;
    }
    // Return solution vector.
    return (Soln); 
}



/**
 * @brief Checks if a value is NAN
 * 
 * @param var Variable we're checking.
 * 
 */
bool NotNAN(double var)
{
    //NaN != NaN.
    if (var != var)
    {
        return (0);
    }
    else
    {
        return (1);
    }
}


/**
 * @brief Returns One if Matrix contains a NaN.
 * 
 * @param mat Matrix we care about. 
 */
bool MatNoNAN(std::vector<std::vector<double>> mat)
{
    for (uint64_t i = 0; i < mat.size(); i++)
    {
        for (uint64_t j = 0; j < mat[0].size(); j++)
        {
            if (!NotNAN(mat[i][j]))
            {
                return (0);
            }
        }
    }
    return (1);
}


/**
 * @brief Returns an identity matrix of size 'size1, size2'. Useful for constructing penalization matrices. 
 * 
 * @param size Function returns a [size1, size2] identity matrix.
 * 
 * @return size1 by size2 identity matrix.
 */
std::vector<std::vector<double>> Eye(uint64_t size1, uint64_t size2)
{
    // Create output, populate with zeros.
    std::vector<std::vector<double>> out(size1, std::vector<double>(size2, 0));

    // Add ones along diagonal.
    for (uint64_t i = 0; i < out.size(); i++)
    {
        for (uint64_t j = 0; j < out[0].size(); j++)
        {
            if (i == j)
            {
                out[i][j] = 1;
            }
        }
    }
    return (out);
}


/**
 * @brief This function adds a wigglyness penalty to a smoothing spline. The derivation
 * for this penalty matrix is given in Wood --- Generalized Additive Models (2006).
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
 */
std::vector<std::vector<double>> AddWigglyPenalty(double lambda, std::vector<std::vector<double>> B)
{
    std::vector<std::vector<double>> out = B;

    // Initialize a vector of zeroes at the correct size (dimensions given in Wood 2006).
    std::vector<std::vector<double>> D = Eye(B.size() - 2, B.size());

    for (uint64_t i = 0; i < D.size(); i++)
    {
        for (uint64_t j = 0; j < D[0].size(); j++)
        {
            if (i == j)
            {
                // Add one along the diagonal.
                D[i][j] = 1;
            }
            if (i + 1 == j)
            {
                // Add -2 on the first off diagonal.
                D[i][j] = -2;
            }
            if (i + 2 == j)
            {
                // Add 1 on the second off diagonal.
                D[i][j] = 1;
            }
        }
    }

    // Matrix multiplication to recover S = DTD.
    D = MatMul(Transpose(D), D);

    // Add penalty to input matrix.
    for (uint64_t i = 0; i < out.size(); i++)
    {
        for (uint64_t j = 0; j < out.size(); j++)
        {
            out[i][j] += pow(lambda, 2) * D[i][j];
        }
    }

    // Returns (Input + lambda^2 DTD).
    return (out);
}



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
 */
bool IsEye(std::vector<std::vector<double>> A, uint64_t l1, uint64_t l2)
{
    for (uint64_t i = 0; i < l1; i++)
    {
        for (uint64_t j = 0; j < l2; j++)
        {
            if (i != j)
            {
                if (A[i][j] != 0)
                { 
                    // If the off-diagonals aren't zero, fail. 
                    return (0);
                }
            }
            if (i == j)
            {
                if (A[i][j] != 1)
                {
                    // If the diagonals aren't one, fail.
                    return (0);
                }
            }
        }
    }
    // If we never failed, pass. 
    return (1);
}


/**
 * @brief Calculates the matrix inverse using Gauss Jordan Elimination. 
 * 
 * @param A Input Matrix.
 */
std::vector<std::vector<double>> Inverse(std::vector<std::vector<double>> A)
{


    if (A.size() != A[0].size())
    {
        std::cout << "error, matrix is not square!"
                  << "\n";
    }

    // Set up matrices to store values //  Identity matrix to augment onto side.
    std::vector<std::vector<double>> IDtoAugment = Eye(A.size(), A[0].size());
    std::vector<std::vector<double>> Inv(A.size(), std::vector<double>((A[0].size() + IDtoAugment[0].size()), 0));
    std::vector<std::vector<double>> C = A;

    // Augment Matrix.
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
    // Augmentation Complete.

    // Define values for the While loop
    uint64_t maxiters = 1000000; // Really big number to avoid a totally fatal infinite loop.
    uint64_t count = 0;          // Count for while loop.

    // Start actual Gaussian elimination process.
    while (count < maxiters)
    {
        for (uint64_t index = 0; index < Inv.size(); index++)
        {
            // Loop over all the indices from zero -> number of rows (not columns; matrix is augmented and no longer square)
            for (uint64_t rowIndex = 0; rowIndex < A.size(); rowIndex++)
            {

                // Divide out to make leading element into one if we're on a pivot.
                if (index == rowIndex)
                {
                    double Divisor = Inv[rowIndex][rowIndex];
                    for (uint64_t k = 0; k < Inv[0].size(); k++)
                    {
                        Inv[rowIndex][k] = Inv[rowIndex][k] / Divisor;
                    }
                }

                // If already zero we skip to avoid pointless work.
                if (Inv[rowIndex][index] != 0) 
                {

                    // Values for computation after quick division by zero check.
                    double currentVal = Inv[rowIndex][index]; 
                    double rowOneVal = Inv[index][index];

                    if (rowOneVal != 0)
                    {
                        if (rowIndex != index)
                        {
                            double correctionFactor = -(currentVal / rowOneVal);

                            // Add a multiple of correctionfactor*rowoneindex + rowindex DO THIS
                            for (uint64_t iter = 0; iter < Inv[0].size(); iter++)
                            {
                                // Adding multiple of one row to another row.
                                Inv[rowIndex][iter] = Inv[rowIndex][iter] + correctionFactor * Inv[index][iter]; // Add row*multiple to row to reduce to zero.
                            }
                        }
                    }
                }
            }
        }

        // Check if LHS is identity. If so, break.
        if (IsEye(Inv, A.size(), A[0].size()) == 1)
        {
            count = maxiters;
        }
        count = count + 1;
    }

    // Push only LHS into a new matrix.
    std::vector<std::vector<double>> Out = A;

    for (uint64_t i = 0; i < Inv.size(); i++)
    {
        for (uint64_t j = 0; j < A[0].size(); j++)
        {
            Out[i][j] = Inv[i][j + A[0].size()];
        }
    }

    // Return LHS as inverse. 
    return (Out);
}

/**
 * @brief Returns the trace for some matrix A.
 * 
 * @param A Input matrix.
 */
double Trace(std::vector<std::vector<double>> A)
{


    double sum = 0;
    for (uint64_t i = 0; i < A.size(); i++)
    {
        for (uint64_t j = 0; j < A[0].size(); j++)
        {
            if (i == j)
            {
                // Sum every diagonal element.
                sum += A[i][j];
            }
        }
    }
    return (sum);
}


/**
 * @brief This returns the trace of the hat matrix in a penalization problem. E.g. for design matrix X, 
 * 
 * X (XTX + lambda^2 D^2)^-1 XT
 * 
 * @param X Design matrix in a regression problem.
 * @param lambda Penalization term. Make zero for a non-penalized problem.
 * 
 */
double TraceHatMatrix(std::vector<std::vector<double>> X, double lambda)
{
    // Set up X, X transpose.
    std::vector<std::vector<double>> XTX = MatMul(Transpose(X), X);

    // Apply Penalty.
    std::vector<std::vector<double>> PenalizedXTX = AddWigglyPenalty(lambda, XTX);

    // Apply final multiplications.
    std::vector<std::vector<double>> InvertedPenXTX = Inverse(PenalizedXTX);
    std::vector<std::vector<double>> TempMat = MatMul(X, InvertedPenXTX);
    std::vector<std::vector<double>> XinvXTXpenX = MatMul(TempMat, Transpose(X));

    // Return.
    return (Trace(XinvXTXpenX));
}

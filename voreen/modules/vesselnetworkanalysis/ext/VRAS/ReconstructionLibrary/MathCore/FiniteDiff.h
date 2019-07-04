
#ifndef __TUBNAVFINITEDIFF_H
#define __TUBNAVFINITEDIFF_H
/*
   The spacing for calculating the first derivative are either uniformly spaced or not
   if they are uniformly spaced then we can use coefficients to calculate the finite differences
   otherwise we should calculate the weights to be used in order to calculate the derivatives

   Calculation of weights in finite difference formulas.
   Based on the approach defined by Bengt Fornberg in
       [1] Generation of Finite Difference Formulas on Arbitrarily Spaced Grids
       [2] Calculation of weights in finite difference formulas.

   The centerline points may not be given with the uniform spacing, therefore the
   finite formulas differences are not 100% applicable, as errors may be introduced. 

   This is a one-dimensional calculation. 
   // ---------------------------------

   d^m / dx^m | x=x0 approx.  sum_{v=0}^n delta^m_ {n, v) f( alpha_v)

   where delta are the weights described by the previous alg. 

   // ------------------------------------- 
   The alpha values are the geodesic distance from the start of the curve 
   
   // The weights is 3-indexed, the first iteration of this will work with an unordered_map of a tuple
   // instead of an array... saves memory, but the constant is slightly most costly. 

*/


#include <type_traits>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <exception>
#include <stdexcept>

#include "../Core/tubNavPoint.h"

namespace tubNav {

   inline namespace FiniteDiffCalculator {



       //**************************************************************************************************
       template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
       static void CalculateWeights(const T x0,const int HighestDerivative /*M*/, const std::vector<T>& alpha /*N & alphas*/, std::vector<T>& fnlWeights ){
            //The pseudo-code is replicated from [1,2]
            //Still need to check the results

            //Only if a safe cast exist, between the desired types
            const T one = static_cast<T>(1);
            const T zero = static_cast<T>(0);

            // The indexing is different, as it is in the pseudocode the
            // range is including the end
            const int M = HighestDerivative ;
            const int N = static_cast<short int>(alpha.size()) - 1;

            T c1 = one; // c1 := 1
            T c4 = alpha.at(0) - x0; // c4 = x(0) - z
            // Initialize weights array to zeros
            T weights[N+1][M+1];
            for(int i = 0; i < N+1; i++)
                for(int j = 0; j <M+1; j++)
                     weights[i][j] = zero;

            weights[0][0] = one;

            for (int i = 1; i <= N; i++){ // for n = 1 to N do
                auto mn = std::min(i, M);
                T c2 = one; // c2 = 1.0
                T c5 = c4;
                c4 = alpha.at(i) - x0;

                for(int j = 0; j <= i -1 ; j++){  // for v := 0 to n-1 do
                     T c3 = alpha.at(i) - alpha.at(j); // c3 := alpha_n - alpha_v
                     c2 = c2 * c3;                     // c2 := c2 * c3

                     if ( j == i-1){ // j. eq. i-1
                         for(int k = mn; k >= 1; k--){
                             weights[i][k] = c1*( k*weights[i-1][k-1] - c5*weights[i-1][k] )/c2;
                         }
                         weights[i][0] = (-c1*c5*weights[i-1][0])/c2;
                     }

                     for(int k = mn; k >= 1; k--){
                         weights[j][k] = (c4*weights[j][k] - k*weights[j][k-1])/c3;
                     }
                     weights[j][0] = c4*weights[j][0] / c3;
                }
                c1 = c2;
            }
            // now save the final weights
            for(int i = 0; i < N + 1; i++){
                fnlWeights.push_back( weights[i][M]);
            }
        }
       // Calculate the alpha based on the complete geodesic distance
       template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
       static void CalculateGeodesicAlpha(std::vector<tubNav::Point<T> >& Points, std::vector<T>& alpha, T (*sqrt)(T)){
           auto distance = static_cast<T>(0);
           alpha.emplace_back(distance);
           for(unsigned int i = 1; i < Points.size(); i++){
              auto c = Points.at(i) - Points.at(i-1);
              distance += sqrt(c.normSquared());
              alpha.emplace_back(distance);
           }
       }
       // Calculate the alpha based only on one dimension of the points
       template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
       static void CalculateGeodesicOneDimensionalAlpha(std::vector<tubNav::Point<T> >& Points, std::vector<T>& alpha, int dim, T (*sqrt)(T)){
           auto distance = static_cast<T>(0);
           alpha.emplace_back(distance);
           for(unsigned int i = 1; i < Points.size(); i++){
              auto c = Points.at(i).get(dim) - Points.at(i-1).get(dim);
              distance += sqrt(c*c);
              alpha.emplace_back(distance);
           }
       }

       // First Derivative from Points at point location i and dimension d
       template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
       static T GetFirstWeightedDerivative(int loc, std::vector<tubNav::Point<T> >& Points, int dim, T (*sqrt)(T)){

           std::vector<double> alphas, weights;
           tubNav::CalculateGeodesicOneDimensionalAlpha(Points, alphas, dim, sqrt);
           tubNav::CalculateWeights( alphas.at(loc) , 1, alphas, weights);

           T result = static_cast<T>(0);
           for(unsigned int i = 0; i < alphas.size(); i++){
               result += weights.at(i)*Points.at(i).get(dim);
           }
           return result;
       }
       //****************************************************************************************************


      const int CentralDifferenceWeights[4][10] = { {-1 ,   0 ,   1, 0 , 0, 0, 0, 0, 0, 2},
                                                    {3 ,  -24,   0, 24,-3, 0, 0, 0 ,0 ,36},
                                                    {-1,   9,   -45, 0, 45,-9, 1, 0, 0,60},
                                                    {105, -1120, 5880, -23520, 0, 23520, -5880, 1120, 105, 29400}
                                             };




       template<typename T, typename Iterator>
       static tubNav::Point<T> GetCentralFiniteDifference( Iterator start, Iterator end, int accuracy ){
           std::size_t size = end - start;
            if ( (size % 2) == 0) {
                std::cout << size << " ..." << size % 2 << std::endl;
                throw std::invalid_argument("Central Finite Difference assumes odd-sized arrays.");
            }
            if ( size < 2 ) {
                throw std::invalid_argument("Minimal central finite difference needs 3 points");
            }
            tubNav::Point<T> point;
            const int idx = (size / 2) -1;

            int i = 0;
            for( ; start != end; ++start, ++i){
                point = point + (*start)*CentralDifferenceWeights[idx][i];
            }
            point = point / CentralDifferenceWeights[idx][9];
            return point;
       }

       // In general, to get the coefficients of the backward approximations, give all odd derivatives listed in the table the
       // opposite sign, whereas for even derivatives the signs stay the same. Given that we are only calculating the
       // first derivative, we just need to do the sign
       const float OneSidedDifferenceWeights[7][10] ={
                                                    {-1 ,  1,     0,    0,      0,     0,     0,   0,  0, 1}, //acc 2 ...idx = 0
                                                    {-3 ,  4,    -1,    0,      0,     0,     0,   0,  0, 2}, //acc 3 .. idx = 1
                                                    {-11,  18,   -9,    2,      0,     0,     0,   0,  0, 6}, //acc 4 .. idx = 2
                                                    {-25,  48,   -36,   16,    -3,     0,     0,   0,  0, 12},//acc 5 .. idx = 3,
                                                    {-137, 300,  -300,  200,   -75,    12,    0,   0,  0, 60},//acc 6 .. idx = 4
                                                    {-147, 360,   450,  400,   -225,   72,   -10,  0,  0, 60}, // acc 7 .. idx = 5
                                                    {-5445,14700,-22050,24500, -18375, 8820, -2450,300,0 ,2100},// acc 8 .. idx = 6
                                                  };

       const float BackWardDifferenceWeights[7][10]  ={{-1 ,  1,     0,    0,      0,     0,     0,   0,  0, 1},//*
                                                     { 1,  -4,     3,    0,      0,     0,     0,   0,  0, 2},//*
                                                     {-2,   9,    -18,  11,      0,     0,     0,   0,  0, 6},//*
                                                     { 3, -16,     36, -48,     25,     0,     0,   0,  0, 12},//*
                                                     {-12, 75,   -200,  300,  -300,    137,    0,   0,  0, 60},//*
                                                     { 10,-72,    225, -400,  -450,   -360,  147,   0,  0, 60},//*
                                                     {-300, 2450, -8820, 18375, -24500, 22050, -14700, 5445,0 , 2100}
                                                   };


       template<typename T, typename Iterator>
       static tubNav::Point<T> GetForwardFiniteDifference( Iterator start, Iterator end, int accuracy ){
           std::size_t size = end - start;
           if ( size < 1 ) throw std::invalid_argument("Forward difference needs at least two points");
           tubNav::Point<T> point;

           const int idx = accuracy -2;

           int i = 0;
           for(; start != end; ++start, ++i){
               point = point + (*start)*OneSidedDifferenceWeights[idx][i];
           }

           point = point / OneSidedDifferenceWeights[idx][9];
           return point;
       }

       template<typename T, typename Iterator>
       static tubNav::Point<T> GetBackwardFiniteDifference( Iterator start, Iterator end, int accuracy ){

           std::size_t size = end - start;
           if ( size < 1 ) throw std::invalid_argument("Backward difference needs at least two points");
           tubNav::Point<T> point;
           const int idx =  accuracy - 2;
           int i = 0;

           for(; start != end; ++start, ++i){
               point = point + (*start)*(BackWardDifferenceWeights[idx][i]);
           }
           point = point / (BackWardDifferenceWeights[idx][9]);

           return point;
       }
  }
}

#endif



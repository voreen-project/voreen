#ifndef __TUBNAV_VECTORMATHOPS_H
#define __TUBNAV_VECTORMATHOPS_H


#include <gsl/gsl_blas.h>
#include "../Core/tubNavPoint.h"
#include <cmath>

namespace tubNav {
   inline namespace MathVectorOps {

        static void SkewSymmetric(const gsl_vector* v, gsl_matrix* Vx){
           gsl_matrix_set_zero(Vx);

           gsl_matrix_set(Vx, 0, 1, -gsl_vector_get(v, 2) );
           gsl_matrix_set(Vx, 1, 0, gsl_vector_get(v, 2) );

           gsl_matrix_set(Vx, 0, 2, gsl_vector_get(v, 1) );
           gsl_matrix_set(Vx, 2, 0, -gsl_vector_get(v, 1) );


           gsl_matrix_set(Vx, 1, 2, -gsl_vector_get(v, 0) );
           gsl_matrix_set(Vx, 2, 1, gsl_vector_get(v, 0) );
        }

        static void RotationAroundAxis( double angle, gsl_vector* axis, gsl_matrix* R){
            // Create the rotation matrix according to Euler-Rodrigues formula.
            double radians = angle * (3.14159/180.0);
            double cosine = cos(radians);
            double sine = sin(radians);
            int degree = 3;

            gsl_matrix_set_zero(R);

            gsl_matrix* I = gsl_matrix_alloc(degree, degree);
            gsl_matrix_set_identity(I);


            gsl_matrix* CrossU = gsl_matrix_alloc(degree, degree);
            gsl_matrix* TensorU = gsl_matrix_alloc(degree, degree);

            SkewSymmetric(axis,  CrossU);

            gsl_matrix_set_zero(TensorU);
            for( int i = 0; i < degree ; i++){
              for(int j = 0; j < degree; j++){
                  gsl_matrix_set(TensorU, i, j, gsl_vector_get(axis, i)* gsl_vector_get(axis,j) );
               }
            }

           gsl_matrix_scale(I, cosine);
           gsl_matrix_scale(CrossU, sine  );
           gsl_matrix_scale(TensorU, (1.0 - cosine));


           gsl_matrix_add(R, I);
           gsl_matrix_add(R, CrossU);
           gsl_matrix_add(R, TensorU);

           gsl_matrix_free(I);
           gsl_matrix_free(CrossU);
           gsl_matrix_free(TensorU);
        }

        inline bool AreEqual(double p1[], double p2[]){
            double error = 0;
            for(int i = 0; i < 3; i++) error += fabs( p1[i] - p2[i]);
            return (error < 0.00001);
        }

        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        static bool IsPointBelow( Point<T> normal, Point<T> pointInPlane, Point<T> pointToCheck){
            /*double n[3] = {normal[0], normal[1], normal[2]};
            double p0[3] = {-pointInPlane[0],-pointInPlane[1],-pointInPlane[2] };
            double p[3] = { pointToCheck[0], pointToCheck[1], pointToCheck[2]};

            double d1  = MathHelper::Dot(p0,n);*/

            Point<T> tmp;
            auto p0 = tmp-pointInPlane;
            double d1  = p0.dot(normal);
            double d2 =  pointToCheck.dot(normal);
            return ( d2 +d1 <= 0);
        }


        static bool AreCloseToEqual(double p1[], double p2[]){
            double error = 0;
            for(int i = 0; i < 3; i++) error += fabs( p1[i] - p2[i]);
            return (error < 0.05);
        }

        static double DotProduct(double p1[], double p2[]){

            Point<double> a(p1[0], p1[1], p1[2], 0);
            Point<double> b(p2[0], p2[1], p2[2], 0);

            return a.dot(b);
        }


        static bool ArePointsCollinear(double x1[], double x2[], double x3[]){

            double triangleArea = x1[0]*( x2[1] - x3[1]) +  x2[0]*(x3[1] - x1[1]) + x3[0]*(x1[1] - x2[1]);
            return ( fabs(triangleArea) < 0.001);
        }

        static double AngleBetweenVectors( double* v1, double*v2, bool degrees){

           gsl_vector *a = gsl_vector_alloc(3);
           gsl_vector *b = gsl_vector_alloc(3);

           double normA = 0;
           double normB = 0;
           for( int i = 0; i < 3; i++){
              gsl_vector_set(a, i, v1[i]);
              normA += v1[i]*v1[i];

              gsl_vector_set(b, i, v2[i]);
              normB += v2[i]*v2[i];
           }

           normA = sqrt(normA);
           normB = sqrt(normB);

           gsl_vector_scale(a, 1.0/normA);
           gsl_vector_scale(b, 1.0/normB);

           double c = 0;
           gsl_blas_ddot(a,b, &c);
           gsl_vector_free(a);
           gsl_vector_free(b);
           // the dot product is 0....
           double origin[3] = {0,0,0};
           double angle =  acos(c);

           if ( ArePointsCollinear(origin, v1,v2)){
               // if the points are collinear, then either they are
               // going in the same direction and the angle is then 0
               // or they are not equal and they are going in different directions....
              if ( AreEqual(v1,v2) || AreCloseToEqual(v1,v2))
                  angle = 0;
              else
                  angle = 3.1415;

           }
           // if the two values are really close, then the result will be 1 and epsilon higher
           // so, if isnan then set to 0 the angle.
           else if ( std::isnan(angle)){
               angle = 0;
           }
           if ( degrees )
               angle = angle * (180.0 / 3.1415);
           return angle;
        }


        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        void MultiplyMatrix(const T* A,const T* B, int n, T* C){
            short int i, j, k;
            T sum = 0;
            for (i = 0; i < n; i++) {
                  for (j = 0; j < n; j++) {
                     sum = 0;
                     for (k = 0; k <  n; k++) {
                        sum = sum + A[i*n + k] * B[k*n + j];
                     }
                     C[i*n + j] = sum;
                  }
            }
        }


        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        void MultiplyMatrixByPoint(T matrix[], T point[], T result[]){

            // point has 4 values
            for(int i =0; i < 4; i++){
                result[i] = 0;
                for(int j = 0; j < 4;j++){
                   result[i] +=  matrix[i*4 + j]*point[j];
                }
            }
        }

        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        void MultiplyMatrixByPoint(T matrix[], Point<T> point, T result[]){

            // point has 4 values
            for(int i =0; i < 4; i++){
                result[i] = 0;
                for(int j = 0; j < 4;j++){
                   if (j < 3)
                      result[i] +=  matrix[i*4 + j]*point[j];
                   else
                      result[i] +=  matrix[i*4 + j];

                }
            }
        }


        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        Point<T> RotatePointAroundAxis( double angle, Point<T> axis, Point<T> point)
        {
            axis.normalize(sqrt);

            gsl_matrix* R = gsl_matrix_alloc(3, 3);
            gsl_vector *a = gsl_vector_alloc(3);
            gsl_vector *p = gsl_vector_alloc(3);

            for(int i = 0; i < 3 ; i++){
               gsl_vector_set(a, i, axis[i]);
               gsl_vector_set(p, i, point[i]);
            }

            RotationAroundAxis( angle, a, R);

            gsl_vector *result = gsl_vector_alloc(3);

            gsl_blas_dgemv(CblasNoTrans, 1.0, R, p, 0, result);

            Point<T> res(gsl_vector_get(result,0), gsl_vector_get(result,1), gsl_vector_get(result,2), 0);

            gsl_matrix_free(R);
            gsl_vector_free(a);
            gsl_vector_free(p);
            gsl_vector_free(result);
            return res;
        }


        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        T AngleBetween3DVectors(Point<T> v1, Point<T> v2){
            double vectorA[3] = {static_cast<double>(v1.getX()), static_cast<double>(v1.getY()), static_cast<double>(v1.getZ()) };
            double vectorB[3] = {static_cast<double>(v2.getX()), static_cast<double>(v2.getY()), static_cast<double>(v2.getZ()) };
            return AngleBetweenVectors(vectorA, vectorB, true);
        }

    }
}



#endif

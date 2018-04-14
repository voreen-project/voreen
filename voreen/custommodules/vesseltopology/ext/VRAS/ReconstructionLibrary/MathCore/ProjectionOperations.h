#ifndef __TUBNAV_PROJECTION_OPER_H
#define __TUBNAV_PROJECTION_OPER_H


#include "../Core/tubNavPoint.h"
#include "../MathCore/VectorOperations.h"

namespace tubNav {

    inline namespace projectionOperations {


        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        void GetCameraMatrix(const Point<T>& position, const Point<T>& focalPoint, const Point<T> viewUp, T cameraMatrix[]){
            // Create Identity

            for(int i = 0; i < 4; i++){
               for(int j = 0; j < 4; j++){
                   if (i == j)
                       cameraMatrix[i*4 +j] =  static_cast<T>(1.0);
                   else
                       cameraMatrix[i*4 + j] = static_cast<T>(0.0);
               }
           }

            Point<T> viewNormal = position - focalPoint;
            viewNormal.normalize(sqrt);
            Point<T> viewSideWays = viewUp.cross(viewNormal);
            viewSideWays.normalize(sqrt);
            Point<T> orthoViewUp = viewNormal.cross(viewSideWays);


            for(int i = 0; i < 3; i++){
                cameraMatrix[0*4 + i] = viewSideWays[i];
                cameraMatrix[1*4 + i] = orthoViewUp[i];
                cameraMatrix[2*4 + i] = viewNormal[i];
            }


            T delta[4], ndelta[4];
            for(int i = 0; i < 3; i++) delta[i] = -position[i];
            delta[3] = 0.0;

            MultiplyMatrixByPoint(cameraMatrix, delta, ndelta);

            cameraMatrix[0*4 + 3] = ndelta[0];
            cameraMatrix[1*4 + 3] = ndelta[1];
            cameraMatrix[2*4 + 3] = ndelta[2];
        }
    }

    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
    void OrthogonalProjectionMatrix(T xmin, T xmax, T ymin, T ymax, T znear, T zfar, T result[]){
        // Set up result as identity
        for(int i = 0; i < 4; i++){
           for(int j = 0; j < 4; j++){
               if (i == j)
                   result[i*4 +j] = 1.0;
               else
                   result[i*4 + j] = 0.0;
           }
       }
       //
       result[0*4 + 0] = 2.0 / (xmax - xmin);
       result[1*4 + 1] = 2.0 / (ymax - ymin);
       result[2*4 + 2] = -2.0 / (zfar - znear);
       result[0*4 + 3] =  -(xmin + xmax)/(xmax - xmin);
       result[1*4 + 3] =  -(ymin + ymax)/(ymax - ymin);
       result[2*4 + 3] =  -(znear + zfar)/(zfar - znear);
    }

    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
    void PointTo2DProjection(T* modelViewMatrix, T* viewport, const Point<T> point, int sizex, int sizey, T* px, T* py){

        T temp[4];
        MultiplyMatrixByPoint(modelViewMatrix, point, temp);
        temp[0] = temp[0] /temp[3];
        temp[1] = temp[1] /temp[3];

        *px =  (temp[0] + 1.0) * (sizex*(viewport[2]-viewport[0])) / 2.0 + sizex*viewport[0];
        *py =   (temp[1] + 1.0) * (sizey*(viewport[3]-viewport[1])) / 2.0 + sizey*viewport[1];
    }

    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
    void GetOrthoMatrix(T orthoMatrix[]){
        // We use the same parameters as VTK by default, this sets
        // the parallel scale to 1.0 and the clipping range from 0.1 to 1000
        T ps = static_cast<T>(1.0);
        T cr[2] = {0.1, 1000.0 };
        OrthogonalProjectionMatrix(-1.0*ps, 1.0*ps, -1*ps, +1*ps, cr[0], cr[1], orthoMatrix);
    }

    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
    Point<T> ProjectPoint(const Point<T>& position, const Point<T>& cameraLoc, const Point<T>& focalPoint, int size[]){
        T modelViewMatrix[16];
        T viewport[4] = {static_cast<T>(0),static_cast<T>(0), static_cast<T>(size[0]), static_cast<T>(size[1]) };
        GetProjectionMatrix(cameraLoc, focalPoint, modelViewMatrix);
        T x,y;
        PointTo2DProjection(modelViewMatrix, viewport, position, size[0], size[1],&x, &y );
        Point<T> projectedPoint(x,y,0, position.getR() );
        return projectedPoint;
    }


    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
    void GetProjectionMatrix(const Point<T> position, const Point<T> focalPoint, T projectionMatrix[]){
        Point<T> tmp = position - focalPoint;
        tmp.normalize(sqrt);

        Point<T> viewUp(0.0, 1.0, 0.0, 0.0);

        if (sinf(acos(tmp[2]) < 0)){
            Point<T> newPoint(0.0, -1.0, 0.0, 0.0);
            viewUp = newPoint;
        }

        T orthoMatrix[16];
        GetOrthoMatrix(orthoMatrix);

        T cameraMatrix[16];
        GetCameraMatrix(position, focalPoint, viewUp, cameraMatrix);
        MultiplyMatrix(orthoMatrix, cameraMatrix, 4, projectionMatrix);

    }

}


#endif

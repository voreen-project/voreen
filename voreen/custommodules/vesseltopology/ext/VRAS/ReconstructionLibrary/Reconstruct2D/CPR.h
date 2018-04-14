#ifndef __TUBNAV__PROJECTEDCPR_H
#define __TUBNAV__PROJECTEDCPR_H

/*
 *   Projected CPR Based on  Kanitsar et al. CPR- Curved Planar Reformation
 *   The input is the segment to project.
 *   and the projection vector.
 *   It basically does a parallel projection of the centerline.
 *   The idea is basically from here to be able to calculate
 *   points to be queried in the volume, and their resulting position.
 *   That way is not dependent on the image type, or how to query it.
 *   Or the interpolation method.
 */


#include "../Core/tubNavSegment.h"
#include "../MathCore/CoordinateFrame.h"
#include "../MathCore/ProjectionOperations.h"
#include "limits.h"
#include "cprinfo.h"

namespace tubNav {
    namespace Reconstruction2D {


        enum class StretchProjectingPlane { X, Y, Z };

        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        CPRInfo CreateProjectedCPR(std::shared_ptr< Segment<T>> segment,const T angle,const T queryStep,const T limitDistance,
                                   const Point<T>& cameraLoc, const Point<T>& focalPoint,int imageSize[]){

             CPRInfo cprInfo;
             cprInfo.stepSize = queryStep;

             std::shared_ptr< CoordinateFrame<T> > frame = std::make_shared<CoordinateFrame<T>>();
             frame->setNavigationCurve(segment);
             frame->InitTraversal();
             double avgDistance = 0;
             while(!frame->IsAtEnd()){
                 std::vector<Point<T>> pointsInAngleDir;
                 // Get the points in the angle direction
                 frame->GetAllPointsInLine(angle, 0, limitDistance, queryStep, pointsInAngleDir);
                 // and in the opposite direction
                 frame->GetAllPointsInLine(angle + 180, queryStep, limitDistance, queryStep, pointsInAngleDir);

                 // for every point we get the projection
                 // It doesn't matter the re-ordering as each point is projected
                 for(auto point: pointsInAngleDir){
                     ReconstructionPixel newReconstructionPixel;
                     newReconstructionPixel.originalLocation = point;
                     newReconstructionPixel.projectedLocation = ProjectPoint(point, cameraLoc, focalPoint, imageSize);
                     cprInfo.cprPoints.push_back(newReconstructionPixel);
                 }
                 auto prev = frame->GetCenterPoint();
                 frame->Next();
                 auto next = frame->GetCenterPoint();
                 avgDistance += next.distance(prev, sqrt);
             }
             cprInfo.avgPointDistance = avgDistance/ frame->GetSize();
             return cprInfo;
         }

        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        CPRInfo CreateStretchedCPR(std::shared_ptr< Segment<T>> segment,const T angle,const T queryStep,const T limitDistance, StretchProjectingPlane projectTo){

            // In Stretched CPR, the projection plane remains static, in order to maintain curvature
            // Therefore we need to find, according to which plane we maintain the curvature , how much to each side needs to be queried

            std::shared_ptr< CoordinateFrame<T> > frame = std::make_shared<CoordinateFrame<T>>();
            frame->setNavigationCurve(segment);
            frame->InitTraversal();

            std::pair<double, double> range( -DBL_MAX, DBL_MAX); // max and min


            while(!frame->IsAtEnd()){

                auto ptOneSide = frame->GetPoint(angle, limitDistance);
                auto ptOtherSide = frame->GetPoint(angle + 180, limitDistance);

                // now check the range according to where it should project to
                // NOTE: code could be simplified by using enum instead of enum class
                // and therefore allowing the cast to int, and use that to get the
                // X,Y,Z using point get(dim)... but this may be more readable


                if ( projectTo == StretchProjectingPlane::X){
                    if (ptOneSide.getX() > range.first ) range.first = ptOneSide.getX();
                    if (ptOneSide.getX() < range.second ) range.second = ptOneSide.getX();

                    if (ptOtherSide.getX() > range.first ) range.first = ptOtherSide.getX();
                    if (ptOtherSide.getX() < range.second ) range.second = ptOtherSide.getX();

                }
                else if ( projectTo == StretchProjectingPlane::Y){
                    if (ptOneSide.getY() > range.first ) range.first = ptOneSide.getY();
                    if (ptOneSide.getY() < range.second ) range.second = ptOneSide.getY();

                    if (ptOtherSide.getY() > range.first ) range.first = ptOtherSide.getY();
                    if (ptOtherSide.getY() < range.second ) range.second = ptOtherSide.getY();
                }
                else if (projectTo == StretchProjectingPlane::Z){
                    if (ptOneSide.getZ() > range.first ) range.first = ptOneSide.getZ();
                    if (ptOneSide.getZ() < range.second ) range.second = ptOneSide.getZ();

                    if (ptOtherSide.getZ() > range.first ) range.first = ptOtherSide.getZ();
                    if (ptOtherSide.getZ() < range.second ) range.second = ptOtherSide.getZ();
                }

                frame->Next();
            }

            // Now we have the range...

            CPRInfo cprInfo;

            frame->InitTraversal();
            int row = 0;
            double avgDistance = 0;
            while(!frame->IsAtEnd()){

               auto center = frame->GetCenterPoint();
               // now we have to define how much to each side... taking into account the
               // distance to each side..
               double distancePositive = limitDistance;
               double distanceNegative = limitDistance;


               // Original positive...
               auto oPositive = frame->GetPoint(angle, limitDistance);
               // Is it closer to max or min ... and then
               // define the distance accordingly
               if ( projectTo == StretchProjectingPlane::X){
                   double distanceToMin = fabs(oPositive.getX() - range.second);
                   double distanceToMax = fabs(oPositive.getX() - range.first);

                   if ( distanceToMin < distanceToMax){ //the positive direction is closer to min
                       distancePositive = distanceToMin + limitDistance;
                       distanceNegative = fabs(range.first - center.getX());
                   }
                   else {
                       distancePositive = distanceToMax + limitDistance;
                       distanceNegative = fabs(range.second - center.getX());
                   }
               }

               if ( projectTo == StretchProjectingPlane::Y){
                   double distanceToMin = fabs(oPositive.getY() - range.second);
                   double distanceToMax = fabs(oPositive.getY() - range.first);

                   if ( distanceToMin < distanceToMax){ //the positive direction is closer to min
                       distancePositive = distanceToMin + limitDistance;
                       distanceNegative = fabs(range.first - center.getY());
                   }
                   else {
                       distancePositive = distanceToMax + limitDistance;
                       distanceNegative = fabs(range.second - center.getY());
                   }
               }

               if ( projectTo == StretchProjectingPlane::Z){
                   double distanceToMin = fabs(oPositive.getZ() - range.second);
                   double distanceToMax = fabs(oPositive.getZ() - range.first);

                   if ( distanceToMin < distanceToMax){ //the positive direction is closer to min
                       distancePositive = distanceToMin + limitDistance;
                       distanceNegative = fabs(range.first - center.getZ());
                   }
                   else {
                       distancePositive = distanceToMax + limitDistance;
                       distanceNegative = fabs(range.second - center.getZ());
                   }
               }





               //********************************************
               std::vector<Point<T>> pointsInAngleDir;
               // Get the points in the angle direction
               frame->GetAllPointsInLine(angle, 0, distancePositive, queryStep, pointsInAngleDir);
               std::vector<Point<T>> pointsInOppositeAngleDir;
               // and in the opposite direction
               frame->GetAllPointsInLine(angle + 180, queryStep, distanceNegative, queryStep, pointsInOppositeAngleDir);
               //Given that the points have to go from one side to the next... we reverse this ...
               std::reverse(pointsInOppositeAngleDir.begin(), pointsInOppositeAngleDir.end());

               std::vector<Point<T>> pointsInLine;
               for(auto point: pointsInOppositeAngleDir)
                   pointsInLine.push_back(point);
               for(auto point: pointsInAngleDir)
                   pointsInLine.push_back(point);

               // Now this has all the points in line... now to define the projection as well...
               int colLoc = 0;

               auto prev = frame->GetCenterPoint();
               frame->Next();
               auto next = frame->GetCenterPoint();

               for(auto point: pointsInLine){
                   ReconstructionPixel newReconstructionPixel;
                   newReconstructionPixel.originalLocation = point;

                   int index = cprInfo.cprPoints.size();

                   std::pair<int,int> key1 = std::make_pair(colLoc, row -1);
                   std::pair<int,int> key2 = std::make_pair(colLoc-1, row);

                   std::pair<int,int> addKey = std::make_pair(colLoc, row);

                   double distance = 0;

                   if ( cprInfo.indices.count(key1) > 0){
                       auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                       distance = prevPoint.projectedLocation.getZ();
                       distance += point.distance(prevPoint.originalLocation, sqrt);
                   }
                   else if (cprInfo.indices.count(key2) > 0 ){
                       // otherwise copy the distance of the one to the side
                       auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                       distance = prevPoint.projectedLocation.getZ();
                   }

                   Point<T> straightCPRPoint(colLoc, row,  distance, point.getR() );

                   newReconstructionPixel.projectedLocation = straightCPRPoint;
                   colLoc++;

                   cprInfo.indices[addKey] = index;
                   cprInfo.cprPoints.push_back(newReconstructionPixel);
               }
               row++;

               avgDistance += next.distance(prev, sqrt);
            }
            cprInfo.avgPointDistance = avgDistance / frame->GetSize();
            return cprInfo;
        }

        // Losses all curvature info.. by setting it straight
        template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
        CPRInfo CreateStraightCPR(std::shared_ptr< Segment<T>> segment,const T angle,const T queryStep,const T limitDistance){

            CPRInfo cprInfo;
            cprInfo.stepSize = queryStep;

            std::shared_ptr< CoordinateFrame<T> > frame = std::make_shared<CoordinateFrame<T>>();
            frame->setNavigationCurve(segment);
            frame->InitTraversal();

            double avgDistance = 0;

            int row = 0;
            while(!frame->IsAtEnd()){
                std::vector<Point<T>> pointsInAngleDir;
                // Get the points in the angle direction
                frame->GetAllPointsInLine(angle, 0, limitDistance, queryStep, pointsInAngleDir);
                std::vector<Point<T>> pointsInOppositeAngleDir;
                // and in the opposite direction
                frame->GetAllPointsInLine(angle + 180, queryStep, limitDistance, queryStep, pointsInOppositeAngleDir);
                //Given that the points have to go from one side to the next... we reverse this ...

                std::reverse(pointsInOppositeAngleDir.begin(), pointsInOppositeAngleDir.end());

                std::vector<Point<T>> pointsInLine;
                for(auto point: pointsInOppositeAngleDir)
                    pointsInLine.push_back(point);
                for(auto point: pointsInAngleDir)
                    pointsInLine.push_back(point);

                // Now this has all the points in line... now to define the projection as well...
                int colLoc = 0;

                auto prev = frame->GetCenterPoint();
                frame->Next();
                auto next = frame->GetCenterPoint();

                for(auto point: pointsInLine){
                    ReconstructionPixel newReconstructionPixel;
                    newReconstructionPixel.originalLocation = point;

                    int index = cprInfo.cprPoints.size();

                    std::pair<int,int> key1 = std::make_pair(colLoc, row -1);
                    std::pair<int,int> key2 = std::make_pair(colLoc-1, row);

                    std::pair<int,int> addKey = std::make_pair(colLoc, row);

                    double distance = 0;

                    if ( cprInfo.indices.count(key1) > 0){
                        auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                        distance = prevPoint.projectedLocation.getZ();
                        distance += point.distance(prevPoint.originalLocation, sqrt);
                    }
                    else if (cprInfo.indices.count(key2) > 0 ){
                        // otherwise copy the distance of the one to the side
                        auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                        distance = prevPoint.projectedLocation.getZ();
                    }

                    Point<T> straightCPRPoint(colLoc, row, distance, point.getR() );
                    newReconstructionPixel.projectedLocation = straightCPRPoint;
                    colLoc++;
                    cprInfo.indices[addKey] = index;
                    cprInfo.cprPoints.push_back(newReconstructionPixel);
                }
                row++;

                avgDistance += next.distance(prev, sqrt);
            }

            cprInfo.avgPointDistance = avgDistance / row;

            return cprInfo;
        }

    }
}

#endif

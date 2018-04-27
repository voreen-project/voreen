#ifndef UNFOLDING_H
#define UNFOLDING_H


//Methdos to generate unfolding of a single segment
//i.e. doesn't take bifurcation into account.

// The unfolding can be done in linear and non-linear rays
// As opposed to the CPR methods, Unfolding methods will rely on
// the queried volume, but also return the queried points so that
// other (registered) volumes can also be queried at the same positions

// In the case of linear and non-linear rays, there are several possible stopping criteria
// & Different possible checks, namely:
// For ray cast:
//     Threshold reached... we define a upper or lower threshold of value, if that value is
//         reached, we stop querying
//     Distance reached...  we define a maximum distance to reach... and only that
//     Combined... if a max distance is reached, and threshold is not obtain, we might define
//                 a default value
//
//
// Values for the ray can also be aggregated ( min, max, avg)  or just the last value maintained
//
// Now the ray can be cast at equal angles... this will not necessarily mean that the arc-distance
// from one point to the next will be the same.. as the cross-section may not be circular...
// therefore we can define an arc distance instead, query at certain angles and if the arc-distance is
// over a certain threshold we need to lower the angle size... ( it is assumed that the angle step is small
// so that the euclidean distance is close to the arc-distance. )
//
// So we have equi-angle and equi-arc types of unfolding
//
// and finally, if we are dealing with linear unfolding, then cross-sections may overlap
// this should also be resolved... which is quite easy by checking the previous planes



#include "../Core/tubNavSegment.h"
#include "../MathCore/CoordinateFrame.h"
#include "../MathCore/ProjectionOperations.h"
#include "../../definitions.h"

#include "limits.h"
#include "cprinfo.h"

namespace tubNav {

namespace  Reconstruction2D {

     enum class AggregationType { LAST, MIN, MAX, AVG };

     // If the threshold type is no, then we don't care about thresholding
     // if its over.. then once we reach over a threshold we stop the aggregation
     // if its under... once the value is under the threshold we stop the aggregation
     // if distance threshold > 0, then that is a maximum distance, if the threshold is not reached
     // then the default value is set in that location....

     enum class ThresholdType { NO, OVERVALUE, UNDERVALUE };

     static std::map<AggregationType, double > Initialization = {{AggregationType::LAST,0.0 },
                                                                 {AggregationType::AVG ,0.0 },
                                                                 {AggregationType::MAX ,-DBL_MAX},
                                                                 {AggregationType::MAX ,DBL_MAX}};

     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     void NextNonContiguousFrame(std::shared_ptr< CoordinateFrame<T> > frame){
         // Check whether there is an overlap
         auto tangent = frame->GetTangentDirection();
         auto center = frame->GetCenterPoint();
         // basically, all points in the next level should be over
         // the coordinate frame found by the next level
         frame->Next();
         while(!frame->IsAtEnd()){
               double nextR = frame->GetCenterPoint().getR();
               int pointsBelow = 0;
               // if points below equals zero, means there won't be any
               // given the radius distance that overlap.. in order
               // to allow for more overlap, this R can either be greater radius by
               // percentage or a fixed value
               double cAngle = 0;
               while( cAngle < 360){
                   auto pt = frame->GetPoint(cAngle, nextR);
                   if ( tubNav::IsPointBelow(tangent, center, pt) ) pointsBelow++;
                   cAngle += 20;
               }

               if ( pointsBelow == 0)
                   break;
               else
                   frame->Next();
         }
     }

     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     itk::ContinuousIndex<double, 3> CheckLocation(Point<T> currentQueryPoint,  ImageType::Pointer volume, AggregationType type, ThresholdType thresholdType, bool* keepGoing,
                                                   T defaultValue, double* aggregatedValue){
         ImageType::PointType pt;
         pt[0] = currentQueryPoint.getX();
         pt[1] = currentQueryPoint.getY();
         pt[2] = currentQueryPoint.getZ();
         itk::ContinuousIndex<double, 3> loc;
         volume->TransformPhysicalPointToContinuousIndex(pt, loc);
         // if its over the distance threshold or outside the volume, what do we do with the
         // aggregated value...
         // max? keep... min? keep ... value threshold not reached? default..
         // avg? keep
         if ( !volume->GetLargestPossibleRegion().IsInside(loc)){
             if ( thresholdType != ThresholdType::NO && type == AggregationType::LAST)
                 *aggregatedValue = defaultValue;
             *keepGoing = false;
         }
         return loc;
     }

     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     void Aggregate(AggregationType type, double value, double* aggregatedValue){
         if ( type == AggregationType::LAST)
             *aggregatedValue = value; // we save the last value we saw in the aggregation
         else if ( type == AggregationType::MAX)
         {
             *aggregatedValue = std::max(value, *aggregatedValue);
         }
         else if ( type == AggregationType::MIN){
             *aggregatedValue = std::min(value, *aggregatedValue);
         }
         else if ( type == AggregationType::AVG){
             *aggregatedValue += value;
         }
     }


     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     void CheckThresholds(T distance, bool* keepGoing, double value, double* aggregatedValue,   AggregationType type, ThresholdType thresholdType, T distanceThreshold = -1, T valueThreshold = -1, T defaultValue = -1){


         if ( thresholdType == ThresholdType::OVERVALUE){
             if ( value > valueThreshold)
                 *keepGoing = false;
         }
         else if ( thresholdType == ThresholdType::UNDERVALUE){
             if ( value < valueThreshold)
                 *keepGoing = false;
         }
         if ( distanceThreshold > 0 && distance > distanceThreshold){
             if ( keepGoing ) // if its over the distance threshold, it would be misleading to send the aggregated value so far
                 *aggregatedValue = defaultValue;
             *keepGoing = false;
         }
     }

     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     CPRInfo EquiAngleLinearUnfolding(std::shared_ptr< Segment<T>> segment, const T angleStep, const T queryStep, ImageType::Pointer volume, AggregationType type, ThresholdType thresholdType,
                                      T distanceThreshold = -1, T valueThreshold = -1, T defaultValue = -1){

         CPRInfo cprInfo;
         cprInfo.stepSize = angleStep;

         std::shared_ptr< CoordinateFrame<T> > frame = std::make_shared<CoordinateFrame<T>>();
         frame->setNavigationCurve(segment);
         frame->InitTraversal();
         int row = 0;
         double avgDistance = 0;
         double totalDistance = 0;

         itk::LinearInterpolateImageFunction<ImageType,double >::Pointer volumeInterpolator =  itk::LinearInterpolateImageFunction<ImageType, double >::New();
         volumeInterpolator->SetInputImage(volume);

         while(!frame->IsAtEnd()){

             T angle = 0;
             int col = 0;
             while( angle < 360){ // This is now the angle that we are querying ..
                 bool keepGoing = true;
                 double endValue = defaultValue;
                 Point<T> lastQueriedPoint;

                 double aggregatedValue = Initialization[type];

                 T distance = 0;
                 int numPoints = 0;

                 while(keepGoing){
                     auto currentQueryPoint = frame->GetPoint(angle, distance);
                     auto loc = CheckLocation(currentQueryPoint, volume, type, thresholdType, &keepGoing, defaultValue, &aggregatedValue);
                     auto value = defaultValue;
                     if (keepGoing)
                         value = volumeInterpolator->EvaluateAtContinuousIndex(loc);

                     distanceThreshold = frame->GetCenterPoint().getR();
                     CheckThresholds(distance, &keepGoing,value, &aggregatedValue, type, thresholdType, distanceThreshold, valueThreshold, defaultValue);

                     if (keepGoing){
                         Aggregate<double>(type, value, &aggregatedValue);
                         lastQueriedPoint = currentQueryPoint;
                         numPoints++;
                     }
                     else { // we don't keep going.. so we can define the end value now...
                         if ( numPoints != 0){
                             // if its avg, we divide by the num of samples
                             if ( type == AggregationType::AVG && defaultValue != aggregatedValue)
                                 endValue = aggregatedValue / static_cast<T>(numPoints);
                             else
                                 endValue = aggregatedValue;
                             // otherwise is already the value that is in aggregated value
                         }
                     }
                     distance += queryStep;
                 }

                 double fDistance = 0;
                 int index = cprInfo.cprPoints.size();
                 std::pair<int,int> key1 = std::make_pair(col, row -1);
                 std::pair<int,int> key2 = std::make_pair(col-1, row);
                 std::pair<int,int> addKey = std::make_pair(col, row);

                 if ( cprInfo.indices.count(key1) > 0){
                     auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                     fDistance = prevPoint.projectedLocation.getZ();
                     fDistance += lastQueriedPoint.distance(prevPoint.originalLocation, sqrt);
                 }
                 else if (cprInfo.indices.count(key2) > 0 ){
                     // otherwise copy the distance of the one to the side
                     auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                     fDistance = prevPoint.projectedLocation.getZ();
                 }


                 Point<T> unfoldingPoint(col, row, avgDistance, frame->GetCenterPoint().distance( lastQueriedPoint, sqrt));
                 cprInfo.AddPixel(lastQueriedPoint, unfoldingPoint,frame->GetCenterPoint(), endValue);
                 cprInfo.indices[addKey] = index;

                 col += 1;
                 angle += angleStep;
             }
             row++;

             // Check whether there is an overlap
             auto center = frame->GetCenterPoint();
             // basically, all points in the next level should be over
             // the coordinate frame found by the next level
             NextNonContiguousFrame(frame);
             auto next = frame->GetCenterPoint();
             avgDistance += next.distance(center, sqrt);
             totalDistance += next.distance(center, sqrt);
         }

         avgDistance /= row;
         cprInfo.avgPointDistance = avgDistance;
         return cprInfo;
    }

     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     T ClosestDistanceToSegment(std::shared_ptr< CoordinateFrame<T> > frame, Point<T> location){
          double minDistance = DBL_MAX;

          frame->InitTraversal();

          while(!frame->IsAtEnd()){

              double distance = location.distance(frame->GetCenterPoint(), sqrt);
              if ( distance < minDistance)
                  minDistance = distance;
              frame->Next();
          }

          return minDistance;
     }

     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     Point<T> GetNextPoint( Point<T> currentPoint, const T sizeStep, std::shared_ptr< CoordinateFrame<T> > frame){
         double dv[3] = {0.0, 0.0, 0.0};
         double SobelKernel[27] = {-1,-2,-1,-2,-4,-2,-1,-2,-1, 0,0,0,0,0,0,0,0,0, 1,2,1,2,4,2,1,2,1};

         for( short int i = -1; i <= 1; i++){
           for( short int j = -1; j <= 1; j++){
             for( short int k = -1; k <= 1; k++){
                int dIndex[3] = { i + 1, j + 1, k + 1};

                Point<T> tmp(currentPoint.getX() + 0.1*i,currentPoint.getY()  +0.1*j, currentPoint.getZ() + 0.1*k, currentPoint.getR());
                T dist = ClosestDistanceToSegment(frame, tmp);

                int positionInImageY = dIndex[2] +  dIndex[0]* 3 + dIndex[1]* 3* 3 ;
                int positionInImageX = dIndex[1] +  dIndex[2]* 3 + dIndex[0]* 3* 3 ;
                int positionInImageZ = dIndex[0] +  dIndex[1]* 3 + dIndex[2]* 3* 3 ;
                dv[0] +=  dist * ( SobelKernel[positionInImageX] * (1.0/27.0));
                dv[1] +=  dist * ( SobelKernel[positionInImageY] * (1.0/27.0));
                dv[2] +=  dist * ( SobelKernel[positionInImageZ] * (1.0/27.0));
             }
           }
         }

         double magnitude = sqrt(dv[0]*dv[0] +  dv[1]*dv[1] + dv[2]*dv[2]);
         dv[0] /= magnitude;
         dv[1] /= magnitude;
         dv[2] /= magnitude;

         Point<T> dir(dv[0], dv[1],dv[2],0);
         auto nextPoint = currentPoint + dir*sizeStep;
         return nextPoint;
     }


     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     CPRInfo EquiAngleNonLinearUnfolding(std::shared_ptr< Segment<T>> segment, const T angleStep, const T queryStep, ImageType::Pointer volume, AggregationType type, ThresholdType thresholdType,
                                      T distanceThreshold = -1, T valueThreshold = -1, T defaultValue = -1){

         CPRInfo cprInfo;
         cprInfo.stepSize = angleStep;
         std::shared_ptr< CoordinateFrame<T> > frame = std::make_shared<CoordinateFrame<T>>();
         std::shared_ptr< CoordinateFrame<T> > secondaryFrame = std::make_shared<CoordinateFrame<T>>();

         frame->setNavigationCurve(segment);
         secondaryFrame->setNavigationCurve(segment);

         frame->InitTraversal();
         int row = 0;

         itk::LinearInterpolateImageFunction<ImageType,double >::Pointer volumeInterpolator =  itk::LinearInterpolateImageFunction<ImageType, double >::New();
         volumeInterpolator->SetInputImage(volume);

         double avgDistance = 0;
         while(!frame->IsAtEnd()){

             T angle = 0;
             int col = 0;
             while( angle < 360){
                 // This is now the angle that we are querying ..
                 bool keepGoing = true;
                 double endValue = defaultValue;
                 Point<T> lastQueriedPoint;

                 double aggregatedValue = Initialization[type];

                 T distance = 0;
                 int numPoints = 0;

                 auto currentQueryPoint = frame->GetPoint(angle, queryStep/2.0); // tiny distance threshold to start going in one direction
                 while(keepGoing){
                     auto loc = CheckLocation(currentQueryPoint, volume, type, thresholdType, &keepGoing, defaultValue, &aggregatedValue);
                     auto value = defaultValue;
                     if (keepGoing)
                         value = volumeInterpolator->EvaluateAtContinuousIndex(loc);

                     distanceThreshold = frame->GetCenterPoint().getR();

                     CheckThresholds(distance, &keepGoing,value, &aggregatedValue, type, thresholdType, distanceThreshold, valueThreshold, defaultValue);

                     if (keepGoing){
                         Aggregate<double>(type, value, &aggregatedValue);
                         lastQueriedPoint = currentQueryPoint;
                         numPoints++;
                     }
                     else { // we don't keep going.. so we can define the end value now...
                         if ( numPoints != 0){
                             // if its avg, we divide by the num of samples
                             if ( type == AggregationType::AVG)
                                 endValue = aggregatedValue / static_cast<T>(numPoints);
                             else
                                 endValue = aggregatedValue;
                             // otherwise is already the value that is in aggregated value
                         }
                     }
                     distance += queryStep;
                     currentQueryPoint = GetNextPoint(currentQueryPoint, queryStep, secondaryFrame );
                 }



                 double fDistance = 0;
                 int index = cprInfo.cprPoints.size();
                 std::pair<int,int> key1 = std::make_pair(col, row -1);
                 std::pair<int,int> key2 = std::make_pair(col-1, row);
                 std::pair<int,int> addKey = std::make_pair(col, row);

                 if ( cprInfo.indices.count(key1) > 0){
                     auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                     fDistance = prevPoint.projectedLocation.getZ();
                     fDistance += lastQueriedPoint.distance(prevPoint.originalLocation, sqrt);
                 }
                 else if (cprInfo.indices.count(key2) > 0 ){
                     // otherwise copy the distance of the one to the side
                     auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                     fDistance = prevPoint.projectedLocation.getZ();
                 }


                 Point<T> unfoldingPoint(col, row,fDistance, frame->GetCenterPoint().distance( lastQueriedPoint, sqrt) );
                 cprInfo.AddPixel(lastQueriedPoint, unfoldingPoint,frame->GetCenterPoint(), endValue);
                 cprInfo.indices[addKey] = index;
                 col += 1;
                 angle += angleStep;
             }
             row++;
             // Check whether there is an overlap
             // basically, all points in the next level should be over
             // the coordinate frame found by the next level
             auto prev = frame->GetCenterPoint();
             frame->Next();
             auto next = frame->GetCenterPoint();
             avgDistance += prev.distance(next, sqrt);
         }
         cprInfo.avgPointDistance = avgDistance / row;
         return cprInfo;
     }




     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     CPRInfo EquiArcLinearUnfolding(std::shared_ptr< Segment<T>> segment, const T arcLength, const T queryStep, ImageType::Pointer volume, AggregationType type, ThresholdType thresholdType,
                                      T distanceThreshold = -1, T valueThreshold = -1, T defaultValue = -1){

         CPRInfo cprInfo;
         cprInfo.stepSize = arcLength;

         std::shared_ptr< CoordinateFrame<T> > frame = std::make_shared<CoordinateFrame<T>>();
         frame->setNavigationCurve(segment);
         frame->InitTraversal();
         int row = 0;
         double avgDistance = 0;
         T epsilon = 0.05* arcLength;
         double angleStep = 3.0;
         double initAngleStep = 3.0;

         itk::LinearInterpolateImageFunction<ImageType,double >::Pointer volumeInterpolator =  itk::LinearInterpolateImageFunction<ImageType, double >::New();
         volumeInterpolator->SetInputImage(volume);

         while(!frame->IsAtEnd()){

             T angle = 0;
             int col = 0;
             Point<T> prevColPoint;

             while( angle < 360){ // This is now the angle that we are querying ..
                 bool keepGoing = true;
                 double endValue = defaultValue;
                 Point<T> lastQueriedPoint;

                 double aggregatedValue = Initialization[type];

                 T distance = 0;
                 int numPoints = 0;

                 while(keepGoing){
                     auto currentQueryPoint = frame->GetPoint(angle, distance);
                     auto loc = CheckLocation(currentQueryPoint, volume, type, thresholdType, &keepGoing, defaultValue, &aggregatedValue);
                     auto value = defaultValue;
                     if (keepGoing)
                         value = volumeInterpolator->EvaluateAtContinuousIndex(loc);

                     distanceThreshold = frame->GetCenterPoint().getR()*1.2;
                     CheckThresholds(distance, &keepGoing,value, &aggregatedValue, type, thresholdType, distanceThreshold, valueThreshold, defaultValue);

                     if (keepGoing){
                         Aggregate<double>(type, value, &aggregatedValue);
                         lastQueriedPoint = currentQueryPoint;
                         numPoints++;
                     }
                     else { // we don't keep going.. so we can define the end value now...
                         if ( numPoints != 0){
                             // if its avg, we divide by the num of samples
                             if ( type == AggregationType::AVG)
                                 endValue = aggregatedValue / static_cast<T>(numPoints);
                             else
                                 endValue = aggregatedValue;
                             // otherwise is already the value that is in aggregated value
                         }
                     }
                     distance += queryStep;
                 }

                 // We only add the point in the case that the distance between continuous points is below a
                 // certain epsilon
                 if ( col == 0){ // we have nothing to compare against, so we only use it

                     double fDistance = 0;
                     int index = cprInfo.cprPoints.size();
                     std::pair<int,int> key1 = std::make_pair(col, row -1);
                     std::pair<int,int> key2 = std::make_pair(col-1, row);
                     std::pair<int,int> addKey = std::make_pair(col, row);

                     if ( cprInfo.indices.count(key1) > 0){
                         auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                         fDistance = prevPoint.projectedLocation.getZ();
                         fDistance += lastQueriedPoint.distance(prevPoint.originalLocation, sqrt);
                     }
                     else if (cprInfo.indices.count(key2) > 0 ){
                         // otherwise copy the distance of the one to the side
                         auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                         fDistance = prevPoint.projectedLocation.getZ();
                     }


                     Point<T> unfoldingPoint(col, row, fDistance, frame->GetCenterPoint().distance( lastQueriedPoint, sqrt) );
                     cprInfo.AddPixel(lastQueriedPoint, unfoldingPoint, frame->GetCenterPoint(), endValue);
                     cprInfo.indices[addKey] = index;

                     col += 1;
                     angle += angleStep;
                     prevColPoint = lastQueriedPoint;
                 }
                 else {
                     double distanceFromPrevCol = prevColPoint.distance(lastQueriedPoint,sqrt);
                     if ( fabs(distanceFromPrevCol - arcLength ) < epsilon){ // if it is bellow the threshold
                         // we add the point...
                         double fDistance = 0;
                         int index = cprInfo.cprPoints.size();
                         std::pair<int,int> key1 = std::make_pair(col, row -1);
                         std::pair<int,int> key2 = std::make_pair(col-1, row);
                         std::pair<int,int> addKey = std::make_pair(col, row);

                         if ( cprInfo.indices.count(key1) > 0){
                             auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                             fDistance = prevPoint.projectedLocation.getZ();
                             fDistance += lastQueriedPoint.distance(prevPoint.originalLocation, sqrt);
                         }
                         else if (cprInfo.indices.count(key2) > 0 ){
                             // otherwise copy the distance of the one to the side
                             auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                             fDistance = prevPoint.projectedLocation.getZ();
                         }




                         Point<T> unfoldingPoint(col, row,fDistance, frame->GetCenterPoint().distance( lastQueriedPoint, sqrt) );
                         cprInfo.AddPixel(lastQueriedPoint, unfoldingPoint, frame->GetCenterPoint(), endValue);
                         cprInfo.indices[addKey] = index;

                         col += 1;
                         angleStep = initAngleStep;
                         angle += angleStep;
                         prevColPoint = lastQueriedPoint;
                     }
                     else {
                          // Otherwise, check if the distance is greater or smaller.
                         // if its greater, we went too far
                         if ( distanceFromPrevCol > arcLength)
                         {
                             if ( angleStep > 0.1)
                                angleStep /= 2.0;
                             angle -= angleStep;
                         }
                         else {
                             if ( angleStep > 0.1)
                                angleStep /= 2.0;
                             angle += angleStep;
                         }
                     }
                 }
             }
             row++;

             // Check whether there is an overlap
             auto center = frame->GetCenterPoint();
             // basically, all points in the next level should be over
             // the coordinate frame found by the next level
             NextNonContiguousFrame(frame);
             auto next = frame->GetCenterPoint();
             avgDistance += next.distance(center, sqrt);
         }

         avgDistance /= row;
         cprInfo.avgPointDistance = avgDistance;
         return cprInfo;
    }


     template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     CPRInfo EquiArcNonLinearUnfolding(std::shared_ptr< Segment<T>> segment, const T arcLength, const T queryStep, ImageType::Pointer volume, AggregationType type, ThresholdType thresholdType,
                                      T distanceThreshold = -1, T valueThreshold = -1, T defaultValue = -1){

         CPRInfo cprInfo;
         cprInfo.stepSize = arcLength;

         std::shared_ptr< CoordinateFrame<T> > frame = std::make_shared<CoordinateFrame<T>>();
         std::shared_ptr< CoordinateFrame<T> > secondaryFrame = std::make_shared<CoordinateFrame<T>>();

         frame->setNavigationCurve(segment);
         secondaryFrame->setNavigationCurve(segment);
         frame->InitTraversal();


         int row = 0;
         double avgDistance = 0;
         T epsilon = 0.05* arcLength;
         double angleStep = 3.0;
         double initAngleStep = 3.0;

         itk::LinearInterpolateImageFunction<ImageType,double >::Pointer volumeInterpolator =  itk::LinearInterpolateImageFunction<ImageType, double >::New();
         volumeInterpolator->SetInputImage(volume);

         while(!frame->IsAtEnd()){

             T angle = 0;
             int col = 0;
             Point<T> prevColPoint;

             while( angle < 360){ // This is now the angle that we are querying ..
                 bool keepGoing = true;
                 double endValue = defaultValue;
                 Point<T> lastQueriedPoint;

                 double aggregatedValue = Initialization[type];

                 T distance = 0;
                 int numPoints = 0;

                 auto currentQueryPoint = frame->GetPoint(angle, queryStep/2.0); // tiny distance threshold to start going in one direction
                 while(keepGoing){
                     auto loc = CheckLocation(currentQueryPoint, volume, type, thresholdType, &keepGoing, defaultValue, &aggregatedValue);
                     auto value = defaultValue;
                     if (keepGoing)
                         value = volumeInterpolator->EvaluateAtContinuousIndex(loc);

                     distanceThreshold = frame->GetCenterPoint().getR();
                     CheckThresholds(distance, &keepGoing,value, &aggregatedValue, type, thresholdType, distanceThreshold, valueThreshold, defaultValue);

                     if (keepGoing){
                         Aggregate<double>(type, value, &aggregatedValue);
                         lastQueriedPoint = currentQueryPoint;
                         numPoints++;
                     }
                     else { // we don't keep going.. so we can define the end value now...
                         if ( numPoints != 0){
                             // if its avg, we divide by the num of samples
                             if ( type == AggregationType::AVG)
                                 endValue = aggregatedValue / static_cast<T>(numPoints);
                             else
                                 endValue = aggregatedValue;
                             // otherwise is already the value that is in aggregated value
                         }
                     }
                     distance += queryStep;
                     currentQueryPoint = GetNextPoint(currentQueryPoint, queryStep, secondaryFrame );
                 }

                 // We only add the point in the case that the distance between continuous points is below a
                 // certain epsilon
                 if ( col == 0){ // we have nothing to compare against, so we only use it

                     double fDistance = 0;
                     int index = cprInfo.cprPoints.size();
                     std::pair<int,int> key1 = std::make_pair(col, row -1);
                     std::pair<int,int> key2 = std::make_pair(col-1, row);
                     std::pair<int,int> addKey = std::make_pair(col, row);

                     if ( cprInfo.indices.count(key1) > 0){
                         auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                         fDistance = prevPoint.projectedLocation.getZ();
                         fDistance += lastQueriedPoint.distance(prevPoint.originalLocation, sqrt);
                     }
                     else if (cprInfo.indices.count(key2) > 0 ){
                         // otherwise copy the distance of the one to the side
                         auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                         fDistance = prevPoint.projectedLocation.getZ();
                     }

                     Point<T> unfoldingPoint(col, row, fDistance, frame->GetCenterPoint().distance( lastQueriedPoint, sqrt) );
                     cprInfo.AddPixel(lastQueriedPoint, unfoldingPoint,frame->GetCenterPoint(), endValue);
                     cprInfo.indices[addKey] = index;

                     col += 1;
                     angle += angleStep;
                     prevColPoint = lastQueriedPoint;
                 }
                 else {
                     double distanceFromPrevCol = prevColPoint.distance(lastQueriedPoint,sqrt);
                     if ( fabs(distanceFromPrevCol - arcLength ) < epsilon){ // if it is bellow the threshold
                         // we add the point...
                         double fDistance = 0;
                         int index = cprInfo.cprPoints.size();
                         std::pair<int,int> key1 = std::make_pair(col, row -1);
                         std::pair<int,int> key2 = std::make_pair(col-1, row);
                         std::pair<int,int> addKey = std::make_pair(col, row);

                         if ( cprInfo.indices.count(key1) > 0){
                             auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key1]];
                             fDistance = prevPoint.projectedLocation.getZ();
                             fDistance += lastQueriedPoint.distance(prevPoint.originalLocation, sqrt);
                         }
                         else if (cprInfo.indices.count(key2) > 0 ){
                             // otherwise copy the distance of the one to the side
                             auto prevPoint = cprInfo.cprPoints[cprInfo.indices[key2]];
                             fDistance = prevPoint.projectedLocation.getZ();
                         }

                         Point<T> unfoldingPoint(col, row, fDistance, frame->GetCenterPoint().distance( lastQueriedPoint, sqrt) );
                         cprInfo.AddPixel(lastQueriedPoint, unfoldingPoint,frame->GetCenterPoint(),  endValue);
                         cprInfo.indices[addKey] = index;

                         col += 1;
                         angleStep = initAngleStep;
                         angle += angleStep;
                         prevColPoint = lastQueriedPoint;
                     }
                     else {
                          // Otherwise, check if the distance is greater or smaller.
                         // if its greater, we went too far
                         if ( distanceFromPrevCol > arcLength)
                         {
                             if ( angleStep > 0.1)
                                angleStep /= 2.0;
                             angle -= angleStep;
                         }
                         else {
                             if ( angleStep > 0.1)
                                angleStep /= 2.0;
                             angle += angleStep;
                         }
                     }
                 }
             }
             row++;
             std::cout << row << "/" << frame->GetSize() << std::endl;

             // Check whether there is an overlap
             auto center = frame->GetCenterPoint();
             // basically, all points in the next level should be over
             // the coordinate frame found by the next level
             frame->Next();
             auto next = frame->GetCenterPoint();
             avgDistance += next.distance(center, sqrt);
         }

         avgDistance /= row;
         cprInfo.avgPointDistance = avgDistance;
         return cprInfo;
    }



}
}

#endif // UNFOLDING_H

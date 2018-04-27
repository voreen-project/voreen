#ifndef __TUBNAVSEGMENT_H
#define __TUBNAVSEGMENT_H

/*
   A segment is defined as the points from the root or
   bifurcation +1 points towards the next bifurcation or end-point.
   
   Local or Global ( segment-wise) information can be accesed from here.
   The tangent, tortuosity, inner & outer angle, as well as linear and
   non-linear ray casting methods are done here. 

   The results may be cached or not. If the representation changes || 
   desired precision || step-size changes then the cached is clear.

   Representation can be set either as Spline or OrderedPointSet.
   The spline can either be a cardinal spline, kochanek-bartels, and finite difference for tangent calculation.
  
   The tangent can be calculated by finite difference,
   by default we use central difference  https://en.wikipedia.org/wiki/Finite_difference_coefficient 
   it assumes uniform spacing. Let's figure this out, as technically the spacing is not uniform when 
   calculating the tangent, therefore there's an error introduced.

   In order to calculate the difference, we use the Bengt Fornberg paper > Calculation of weights in finite difference
   formula.

   Then we have the ray cast, We use Wang et al. approach for the computation of rotation minimizing frames...
   Which is more stable than Frechet or Klok. 
   Want approach needs an initial right-handed orthonormal frame... The calculation of this will be detailed later.
   It needs to be reproducible, so a random vector is not possible. And pre-selecting a single vector may cause errors..
   

   We also calculate Tortuosity in the segments with the measures used
   "Measuring Tortuosity of the Intracerebral Vasculature from MRA Images" Bullitt et al. 

   Let's also define Frechet Distance Variants, especially the Euclidean Transforms Invariants.    

   The geodesic distance as well. 
*/


#include "tubNavPoint.h"
#include "../MathCore/FiniteDiff.h"

#include <vector>
#include <iterator>
#include <memory>
#include <type_traits>
#include <bitset>
#include <cmath>

namespace tubNav {


   template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
   class SegmentIterator ;

   template<typename T,
            typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
   class Segment {

          friend class SegmentIterator<T>;
      public:

          Segment(): _fixedLocationIndex(0),_evaluationLocation(0),_stepSize(0.02), _typeCurve(Representation::POINTS),_modified(false),_accuracy(3),
                     _useUniformSpacing(true)
          {

          }


          void SetToStart(){
              _fixedLocationIndex = 0;
              _evaluationLocation = 0;
          }

          void AddPoint(Point<T> point){
              _pointsInSegment.emplace_back(point);
              _insideOthers.push_back(false);
              _modified = true;
          }

          void AddPoints(const std::vector<Point<T>>& points){
              _pointsInSegment = points;
              _insideOthers.clear();
              for(unsigned int i =0; i < points.size(); i++)
                  _insideOthers.push_back(false);
              _modified = true;
          }


          void GetAllPoints(std::vector<Point<T>>& points){
              for(unsigned int i =0; i < _pointsInSegment.size(); i++){
                  points.push_back(_pointsInSegment.at(i));
              }
          }

          void Append(const Segment& o){


          }
          Point<T> GetCurrentPoint(){
              Point<T> point;
              if ( _typeCurve == Representation::SPLINE)
                  point = GetSplinePoint(_evaluationLocation);
              if ( _typeCurve == Representation::POINTS)
                  point = GetVectorPoint(_fixedLocationIndex);
              return point;
          }

          bool GetCurrentPointInsideProperty(){
              return _insideOthers.at(_fixedLocationIndex);
          }

          void SetCurrentPointInsideProperty(bool isInside){
               _insideOthers[_fixedLocationIndex] = isInside;
          }

          Point<T> GetNextPoint(){
              Point<T> point;
              if ( _typeCurve == Representation::SPLINE)
                  point = GetSplinePoint(_evaluationLocation + _stepSize);
              if ( _typeCurve == Representation::POINTS)
                  point = GetVectorPoint(_fixedLocationIndex + 1);
              return point;
          }

          Point<T> GetPrevPoint(){
              Point<T> point;
              if ( _typeCurve == Representation::SPLINE)
                  point = GetSplinePoint(_evaluationLocation - _stepSize);
              if ( _typeCurve == Representation::POINTS)
                  point = GetVectorPoint(_fixedLocationIndex - 1);
              return point;
          }

          Point<T> GetTangentAtCurrentLocation(){

              Point<T> tangent;
              if ( _useUniformSpacing){
                  tangent = GetUniformSpacingTangent();
              }
              else {
                  tangent = GetArbitrarySpacingTangent();
              }
              return tangent;
          }

          void SetRepresentationToSpline(){
              _typeCurve = Representation::SPLINE;
              _modified = true;

          }
          void SetRepresentationToPoints(){
              _typeCurve = Representation::POINTS;
              _modified = true;
          }

          std::size_t GetSize(){
              std::size_t size = 0;
              if ( _typeCurve == Representation::SPLINE && !_pointsInSegment.empty()){
                   size = 1;
              }
              if (_typeCurve == Representation::POINTS)
                  size = _pointsInSegment.size();
              return size;
          }

          void UpdateInfo(){
              // should update info
          }

          void Modified(bool modified){ _modified = modified;}
          bool IsModified(){ return _modified; }

          Point<T> GetStart(){
             Point<T> point;
              if (_typeCurve == Representation::POINTS)
                  point = GetVectorPoint(0);
              if (_typeCurve == Representation::SPLINE)
                  point = GetSplinePoint(0);
              return point;
          }

          Point<T> GetEnd(){
              Point<T> point;

              if (_typeCurve == Representation::POINTS){
                  point = GetVectorPoint(GetSize()-1);
              }
              if (_typeCurve == Representation::SPLINE)
                  point = GetSplinePoint(GetSize());
              return point;
          }

          void MoveForward(){
              if (_typeCurve == Representation::POINTS ){
                  int n =  _pointsInSegment.size();
                  if ( _fixedLocationIndex + 1 <= n)
                      _fixedLocationIndex += 1;
              }
              if (_typeCurve == Representation::SPLINE){
                  if ( _evaluationLocation + _stepSize < 1.0) _evaluationLocation += _stepSize;
                  else _evaluationLocation = 1.0;
              }
          }

          enum class Representation { SPLINE, POINTS};

       private:

           Point<T> GetUniformSpacingPointTangent(int loc){
               bool success = true;

               Point<T> tangent = GetUniformSpacingPointTangentWithAccuracy(loc, _accuracy, &success );

               if (!success){
                   tangent = GetUniformSpacingPointTangentWithAccuracy(loc, 2, &success );
               }
               return tangent;
           }

           Point<T> GetUniformSpacingPointTangentWithAccuracy(int loc, int accuracy, bool* success){
               Point<T> tangent;
               int totalPointstoSide = accuracy;
               if ( accuracy % 2 == 0) totalPointstoSide -= 1; // The central differences uses 3,5,7, and 9
               totalPointstoSide /= 2;
               if (totalPointstoSide == 0) totalPointstoSide+=1;
               // Check if it can be calculated
               int n = _pointsInSegment.size();

               if ( loc - totalPointstoSide >= 0 && loc + totalPointstoSide +1 < n ){
                       typename std::vector<Point<T> >::iterator start = _pointsInSegment.begin() + ( loc - totalPointstoSide);
                       typename std::vector<Point<T> >::iterator end = _pointsInSegment.begin() + ( loc  + totalPointstoSide +1);
                       tangent = GetCentralFiniteDifference<T>(start, end, accuracy);
               }
               else if ( loc - accuracy < 0 && loc + accuracy -1 < n && loc >= 0){ // If we don't have enough points to the left, then we do forward differencing
                       typename std::vector<Point<T> >::iterator start = _pointsInSegment.begin() + ( loc );
                       typename std::vector<Point<T> >::iterator end = _pointsInSegment.begin() + (loc + accuracy);
                       tangent = GetForwardFiniteDifference<T>(start, end, accuracy);
               }
               else if ( loc + accuracy >= n && loc - accuracy >= 0 && loc < n){ // If we don't have enough points to the right, then we do backward differencing

                       typename std::vector<Point<T> >::iterator start = _pointsInSegment.begin() + (loc - accuracy);
                       typename std::vector<Point<T> >::iterator end = _pointsInSegment.begin() + ( loc );
                       tangent = GetBackwardFiniteDifference<T>(start, end, accuracy);
               }
               else {
                       if ( loc == 1){
                             tangent = _pointsInSegment.at(1) - _pointsInSegment.at(0);
                       }
                  *success = false;
               }
               return tangent;
           }

           Point<T> GetUniformSpacingTangent(){
                  Point<T> tangent;
                  if (_typeCurve == Representation::POINTS ){
                      tangent = GetUniformSpacingPointTangent(_fixedLocationIndex);
                  }

                  if (_typeCurve == Representation::SPLINE){

                  }
                  return tangent;
           }

          Point<T> GetNextTangent(){
               Point<T> tangent;
               if (_useUniformSpacing){
                   tangent = GetUniformSpacingPointTangent(_fixedLocationIndex + 1);
               }
               return tangent;
          }

           Point<T> GetArbitrarySpacingTangent(){
               Point<T> tangent;

               return tangent;
           }

           Point<T> GetVectorPoint(const int idx) const{
              if (idx >= static_cast<int>(_pointsInSegment.size())){
                  std::cout << "Getting last point" << std::endl;
                  return _pointsInSegment.back();
              }
              return _pointsInSegment.at(idx);
           }

           Point<T> GetSplinePoint(const double evaluateAt ) const {
               Point<T> point;
               return point;
           }

           void MoveForward(const int movement){
               if (_typeCurve == Representation::POINTS ){
                   int n =  _pointsInSegment.size();
                   if ( _fixedLocationIndex + movement <= n)
                       _fixedLocationIndex += movement;
               }
               if (_typeCurve == Representation::SPLINE){
                   if ( _evaluationLocation + movement*_stepSize < 1.0) _evaluationLocation += movement*_stepSize;
                   else _evaluationLocation = 1.0;
               }
           }

           void MoveBackward(){
               if (_typeCurve == Representation::POINTS ){
                   if ( _fixedLocationIndex - 1 >= 0)
                       _fixedLocationIndex -= 1;
               }
               if (_typeCurve == Representation::SPLINE){
                   if ( _evaluationLocation - _stepSize >= 0 ) _evaluationLocation -= _stepSize;
                   else _evaluationLocation = 0.0;
               }
           }


           void MoveBackward(const int movement){
               if (_typeCurve == Representation::POINTS ){
                   if ( _fixedLocationIndex - movement >= 0)
                       _fixedLocationIndex -= movement;
               }
               if (_typeCurve == Representation::SPLINE){
                   if ( _evaluationLocation - movement*_stepSize >= 0 ) _evaluationLocation -= movement*_stepSize;
                   else _evaluationLocation = 0.0;
               }
           }
           void ReturnToStart(){
               if (_typeCurve == Representation::POINTS ) _fixedLocationIndex = 0;
               if (_typeCurve == Representation::SPLINE ) _evaluationLocation = 0;
           }

           bool IsAtEnd(){
               bool atEnd = false;
               if (_typeCurve == Representation::POINTS ) atEnd = (GetCurrentPoint() == GetEnd());
               if (_typeCurve == Representation::SPLINE ) atEnd = fabs(_evaluationLocation - 1.0) < std::numeric_limits<T>::epsilon();
               return atEnd;
           }

           std::vector<Point<T>> _pointsInSegment;
           std::vector<bool> _insideOthers;

           // If the representation is an ordered point set, then we have a fixed location
           int _fixedLocationIndex;
           // If the representation is a spline then evaluate at u
           double _evaluationLocation; // from 0 ... 1
           double _stepSize;
           Representation _typeCurve;
           bool _modified;
           int  _accuracy; // Accuracy for the tangent function, i.e. how many points should be used surrounding it
           bool _useUniformSpacing;

   };

   template<typename T, typename U >
   class SegmentIterator : public std::iterator<std::random_access_iterator_tag, Segment<T> >{


      public:
          SegmentIterator(std::shared_ptr<Segment<T> > segment) :_segment(segment) { }
          void SetCurve(std::shared_ptr<Segment<T> > segment){ _segment = segment; InitTraversal(); }
          SegmentIterator& operator++(){ _segment->MoveForward(); return *this; } // prefix
          SegmentIterator& operator--(){ _segment->MoveBackward(); return *this; }
          SegmentIterator& operator+=(const int step){ _segment->MoveForward(step); return *this;}
          SegmentIterator& operator-=(const int step){ _segment->MoveBackward(step); return *this;}


          void InitTraversal(){ _segment->ReturnToStart(); }
          bool IsAtEnd(){ return _segment->IsAtEnd(); }
          bool IsAtStart(){ return _segment->GetCurrentPoint() == _segment->GetStart();  }
          bool IsInside(){ return _segment->GetCurrentPointInsideProperty(); }

          bool operator==(const SegmentIterator& rhs) const { return rhs.GetCurrentPoint() == _segment->GetCurrentPoint();   }

          std::size_t GetSize(){ return _segment->GetSize(); }
          Point<T> GetTangent(){ return _segment->GetTangentAtCurrentLocation(); }
          Point<T> GetTangentDirection(){
                Point<T> tangent = _segment->GetTangentAtCurrentLocation();
                // now we normalize the tangent
                tangent.normalize(sqrt);
                return tangent;
          }


          bool IsNextEnd(){
              return _segment->GetNextPoint() == _segment->GetEnd();
          }
          bool IsPrevStart(){ return _segment->GetPrevPoint() == _segment->GetStart();  }
          Point<T> GetPoint(){ return _segment->GetCurrentPoint(); }
          Point<T> GetNextPoint(){ return _segment->GetNextPoint(); }
          Point<T> GetNextTangent(){
              return _segment->GetNextTangent();
          }

          Point<T> GetNextTangentDirection(){
              Point<T> tangent = _segment->GetNextTangent();
              if ( !tangent.IsZeroth())
                 tangent.normalize(sqrt);

              return tangent;
          }

          Point<T> GetPrevPoint(){ return _segment->GetPrevPoint(); }
      private:
          std::shared_ptr<Segment<T> > _segment;
   };

}
#endif

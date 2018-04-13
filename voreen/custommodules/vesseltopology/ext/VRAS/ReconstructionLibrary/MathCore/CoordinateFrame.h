#ifndef __TUBNAV_COORDINATE_FRAME_H
#define __TUBNAV_COORDINATE_FRAME_H


/* The coordinate frame used for the navigation of a segment.
   Based on the paper "Computation of Rotation Minimizing Frames" by Wang et. al
   The original frenet frame defined by frenet-serra formulas are undefined on areas where no curvature exists. In order to avoid this
   an extended frenet-frame as proposed by Klok et al. was initially used, this caused for consecutive frames
   to not rotate too much and also defined in areas where the curvature was not existent.
   Even so, the rotation along the centerline on points "far" apart was to great.
   The current approach fixes those issues, but has one issue. It depends on an
   initial frame, if this frame doesn't exist, then we define our own.

   A random vector is not appropiate as it is not reproducible.
   We need a vector that is in the first coordinate frame, for that we do
   the cross-product of the tangent with a specific vector. This vector
   is based on the coordinate axes. It may happen that the tangent and the
   chosen axis is parallel, and therefore the zero vector is generated. If this happens
   then a next axis is chosen. A default vector can be chosen as well, as it the current segment
   may come from a previous one, and the minimal rotation may cause an offset.
*/

#include "../Core/tubNavSegment.h"

#include <memory>
#include <utility>
#include <vector>
#include <type_traits>

#ifndef __TUBNAV_VECTORMATHOPS_H
    #include "../MathCore/VectorOperations.h"
#endif

namespace tubNav {

     template<typename T,
         typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type >
     class CoordinateFrame {
      public:

         enum class DefaultVectorDirection { XYZ, XZY, YXZ, YZX, ZXY, ZYX};

         CoordinateFrame(): curve(nullptr), dir(DefaultVectorDirection::XYZ ){ }

         CoordinateFrame(DefaultVectorDirection direction): curve(nullptr), dir(direction){ }

         CoordinateFrame(std::shared_ptr<Segment<T> >  segment, DefaultVectorDirection direction =DefaultVectorDirection::XYZ ): dir(direction){
            setNavigationCurve(std::move(segment));
         }

         void SetReferenceVector( Point<T> u0 ){
             _u0 = u0;
         }

         std::shared_ptr<Segment<T> > GetNavigationCurve(){
             return curve;
         }


         void setNavigationCurve(std::shared_ptr<Segment<T> > segment, bool debug = false){
             // We need a copy to mantain here and apply the operations
             curve = segment; //std::move(segment);
             CalculateAllDoubleReflectionReferenceVectors(debug);
         }

         int GetSize(){
             return doubleReflectionReferenceVectors.size();
         }

         Point<T> GetCenterPoint() const {
             SegmentIterator<T> it(curve);
             it.InitTraversal();
             it += _idx;
             return it.GetPoint();
         }

         Point<T> GetTangent() const {
             SegmentIterator<T> it(curve);
             it.InitTraversal();
             it += _idx;
             return it.GetTangent();
         }

         Point<T> GetTangentDirection() const {
             SegmentIterator<T> it(curve);
             it.InitTraversal();
             it += _idx;
             return it.GetTangentDirection();
         }

         Point<T> GetNextCenterPoint() const {
             SegmentIterator<T> it(curve);
             it.InitTraversal();
             it += (_idx+1);
             return it.GetPoint();
         }

         Point<T> GetPoint(const T degrees,const T distance) const {
             SegmentIterator<T> it(curve);
             it.InitTraversal();
             it += _idx;
             Point<T> p = doubleReflectionReferenceVectors.at(_idx)*distance;
             Point<T> rotatedDir = MathVectorOps::RotatePointAroundAxis(degrees, it.GetTangentDirection(), p);
             Point<T> result = it.GetPoint() + rotatedDir;
             return result;
         }


         Point<T> GetNextPoint(const T degrees,const T distance) const {
             SegmentIterator<T> it(curve);
             it.InitTraversal();
             it += (_idx+1);
             Point<T> p = doubleReflectionReferenceVectors.at(_idx+1)*distance;
             Point<T> rotatedDir = MathVectorOps::RotatePointAroundAxis(degrees, it.GetTangentDirection(), p);
             Point<T> result = it.GetPoint() + rotatedDir;
             return result;
         }

         void GetAllPointsInLine(const T degrees, const T startDistance,  const T limitDistance, const T stepSize, std::vector<Point<T>>& allPointsInLine){
             SegmentIterator<T> it(curve);
             it.InitTraversal();
             it += _idx;
             Point<T> p = doubleReflectionReferenceVectors.at(_idx)*1.0;
             Point<T> rotatedDir = RotatePointAroundAxis(degrees, it.GetTangentDirection(), p);

             T distance = startDistance;
             while( distance <= limitDistance){
                 Point<T> newPoint = it.GetPoint() + rotatedDir*distance;
                 distance += stepSize;
                 allPointsInLine.push_back(newPoint);
             }
         }

         Point<T> GetPoint(const int loc, T degrees, T distance){
             SegmentIterator<T> it(curve);
             it += loc;
             Point<T> p = doubleReflectionReferenceVectors.at(loc)*distance;
             Point<T> rotatedDir = RotatePointAroundAxis(degrees, it.GetTangentDirection(), p);
             Point<T> result = it.GetPoint() + rotatedDir;
             return result;
         }

         Point<T> GetLastReferenceVector(){
             Point<T> last = doubleReflectionReferenceVectors.at(doubleReflectionReferenceVectors.size() -1);
             return last;
         }

         void InitTraversal(){ _idx = 0;}
         void SetToEnd(){ _idx = GetSize() -1; }
         void Next(){ _idx++;}
         void Back(){ _idx--;}
         bool IsAtStart(){ return _idx == 0; }
         bool IsAtEnd(){  return _idx == GetSize(); }
         bool IsNextEnd(){  return _idx == GetSize()-1; }

       private:
            void GenerateInitialReferenceVector(){
                // The initial vector, depends on the initial vector direction
                SegmentIterator<T> it(curve);
                it.InitTraversal();

                auto initialTangent = it.GetTangentDirection();
                const int order = static_cast<int>(dir);
                Point<T> possibleReferenceVector;
                for(int i = 0; i < 3; i++){
                    Point<T> testVector =  axes[orderAxes[order][i]];
                    possibleReferenceVector = testVector.cross(initialTangent);
                    if ( possibleReferenceVector.normSquared() > std::numeric_limits<T>::epsilon()*10 ){
                        break;
                    }
                }
                _u0 = possibleReferenceVector;
            }


            //enum class DefaultVectorDirection { XYZ, XZY, YXZ, YZX, ZXY, ZYX}
            const double axes[3][3]   =  { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} };
            const int orderAxes[6][3] =  { {0,1,2}, {0,2,1}, {1,0,2}, {1, 2,0}, {2,0,1}, {2,1,0} };

            void CalculateAllDoubleReflectionReferenceVectors(bool debug=false){

                 if ( debug ){
                    std::cout << "Is u0? " << _u0.IsZeroth() << std::endl;
                 }

                 if ( _u0.IsZeroth() ) // The initial u0 is not fixed, so we calculate it according to the default vector directions
                 {
                    if ( debug) std::cout << "Generating initial vector? " << std::endl;
                    GenerateInitialReferenceVector();
                 }

                 _idx = 0;
                 doubleReflectionReferenceVectors.clear();
                 doubleReflectionReferenceVectors.emplace_back(_u0);

                 SegmentIterator<T> it(curve);
                 if ( debug) std::cout << "Size of curve? " << curve->GetSize() << std::endl;

                 it.InitTraversal();
                 int i = 0;
                 while(!it.IsAtEnd()){ // Just before the end, basically n-1

                     auto v1 = it.GetNextPoint() - it.GetPoint();
                     v1.normalize(sqrt);
                     auto c1 = v1.dot(v1);
                     auto r_i = doubleReflectionReferenceVectors.at(i);
                     auto v1dotRi = v1.dot(r_i);

                     auto r_iL = r_i -   v1*(2.0/c1)*v1dotRi;

                     auto t_i = it.GetTangentDirection();
                     auto v1dotTi = v1.dot(t_i);

                     auto t_iL = t_i -  v1*(2.0/c1)*v1dotTi;

                     auto t_i1 = it.GetNextTangentDirection();

                     auto v2 = t_i1 - t_iL;  // v2 = t_{i+1} - t_i^l
                     v2.normalize(sqrt);
                     auto c2 =  v2.dot(v2);
                     auto v2DotR_iL = v2.dot(r_iL);
                     auto r_i1 = r_iL - v2*(2.0/c2)*v2DotR_iL;
                     r_i1.normalize(sqrt);
                     doubleReflectionReferenceVectors.emplace_back(r_i1);
                     ++it;
                     ++i;
                 }
            }

           // Keep ownership of a shared pointer
           std::shared_ptr<Segment<T> > curve;
           std::vector<Point<T> > doubleReflectionReferenceVectors;
           Point<T> _u0;
           int _idx; // location inside curve
           DefaultVectorDirection dir;
     };
}

#endif

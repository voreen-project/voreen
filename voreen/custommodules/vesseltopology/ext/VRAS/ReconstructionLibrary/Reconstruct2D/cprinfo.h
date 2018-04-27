#ifndef CPRINFO_H
#define CPRINFO_H

// The CPR Info, its the wrapper
// That contains the result from the Projected CPR
// The straight CPR, and the Stretched CPR
// Basically, the points ( with their physical location in the world)
//
#include "../Core/tubNavPoint.h"
#include <vector>
#include <map>

namespace tubNav {
    namespace Reconstruction2D {

       struct ReconstructionPixel{
           Point<double> originalLocation;
           Point<double> centerlineOrigin;
           Point<double> projectedLocation;
           double value;
       };

       class CPRInfo {

            public:
               std::vector<ReconstructionPixel> cprPoints;
               double stepSize; // step size in ray... columns
               double avgPointDistance;
               int dimX;
               int dimY;
               void AddPixel(Point<double> original, Point<double> projected, Point<double> centerlineOrigin, double value){
                   ReconstructionPixel newReconstructionPixel;
                   newReconstructionPixel.originalLocation = original;
                   newReconstructionPixel.projectedLocation = projected;
                   newReconstructionPixel.centerlineOrigin = centerlineOrigin;
                   newReconstructionPixel.value = value;
                   cprPoints.push_back(newReconstructionPixel);

               }
               //TODO- missing clean method for CPR Info

               // When we push, we also write which index is the last value
               // allowing us to find it in the vector & look at the distance ... b/w rows
               // in the case of projected cpr, more than one point should be set,
               // but projected doesn't need that, so only for the others
               std::map< std::pair<int,int>, int> indices;
       };
    }
}
#endif // CPRINFO_H
